"""Independent archive parity against compact real-data fixtures, no API internals."""
from pathlib import Path
from datetime import datetime,timezone
import argparse,json,sys,hashlib
import numpy as np,pandas as pd

R=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser();ap.add_argument('--software-root',type=Path,default=R);ap.add_argument('--fixtures-dir',type=Path,default=R/'validation/fixtures');ap.add_argument('--out',type=Path,default=R/'validation/replay_results.json');args=ap.parse_args()
sys.path.insert(0,str(args.software_root))
from bulk_score_tools import analyze_bulk_scores
implementation_hashes_before={str(p.relative_to(args.software_root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (args.software_root/'bulk_score_tools').rglob('*.py')}
F=args.fixtures_dir
expected=pd.read_csv(F/'expected_correlations.csv');sig=json.loads((F/'fixed_signatures_P1_P6.json').read_text(encoding='utf-8'))
maxgap={'raw':0.,'known_loading':0.,'joint_ols':0.,'archived_raw':0.,'archived_known_loading':0.,'archived_joint_ols':0.,'gene_lambda':0.,'weighted_index':0.,'conditional_support':0.,'unit_or_order_invariance':0.}
counts={'archive_comparable_calls':0,'strict_zero_failures':0,'explicit_complete_case_calls':0,'six_pair_leakage_calls':0,'gene_lambda_rows_checked':0}
diff_records=[];context={}
def check(a,b,key,tol=2e-11):
 if a is None:a=np.nan
 a=np.asarray(a,dtype=float);b=np.asarray(b,dtype=float)
 for i,(av,bv) in enumerate(zip(a.reshape(-1),b.reshape(-1))):
  diff_records.append({**context,'check_type':key,'component_index':i,'observed':float(av),'expected':float(bv),'absolute_difference':float(abs(av-bv)) if np.isfinite(av) and np.isfinite(bv) else np.nan,'tolerance':tol})
 np.testing.assert_allclose(a,b,atol=tol,rtol=0,equal_nan=True)
 finite=np.isfinite(a)&np.isfinite(b)
 if finite.any():maxgap[key]=max(maxgap[key],float(np.max(abs(a-b)[finite])))
def direct_r(a,b):
 if len(a)<3 or min(np.ptp(a),np.ptp(b))<=1e-12:return np.nan
 return float(np.corrcoef(a,b)[0,1])
def compute_expected_on_mask(bulk,q,pair,alpha,method,mask):
 s=[];logs=np.log2(q.loc[mask,['Fibroblast','Macrophage']].to_numpy())
 for side in ['F','M']:
  sd=sig[pair][side];v=np.log2(1e6*bulk.loc[mask,sd['genes']].to_numpy()+alpha)
  s.append(v@np.asarray(sd['weights'])+sd['offset'])
 s=np.column_stack(s)
 if method=='known_loading':corrected=s-logs*np.asarray([sum(sig[pair][side]['weights']) for side in ['F','M']])
 else:
  z=np.column_stack((np.ones(len(s)),logs));corrected=s-z@np.linalg.lstsq(z,s,rcond=None)[0]
 return direct_r(s[:,0],s[:,1]),direct_r(corrected[:,0],corrected[:,1])

for case,g in expected.groupby('case',sort=False):
 bulk=pd.read_csv(F/case/'bulk_probability.csv',index_col=0,float_precision='round_trip')
 q=pd.read_csv(F/case/'rna_q.csv',index_col=0,float_precision='round_trip')
 assert len(bulk)==len(q)==118
 if case.startswith('Simkin'):assert q.columns.tolist()==['Fibroblast','Macrophage','Malignant']
 else:assert q.columns.tolist()==['Fibroblast','Macrophage','Epithelial']
 for (pair,alpha),gg in g.groupby(['pair_id','alpha_CPM'],sort=False):
  eraw=gg[gg.method=='raw'].iloc[0]
  for method in ['known_loading','joint_ols']:
   context={'case':case,'pair_id':pair,'alpha':float(alpha),'method':method}
   ecorr=gg[gg.method==method].iloc[0]
   call=dict(signatures=pair,method=method,bulk_scale='probability',alpha=float(alpha),q_provenance='known',compute_leakage=False)
   if not eraw.full_bulk_score_support:
    try:analyze_bulk_scores(bulk,q,**call)
    except ValueError:counts['strict_zero_failures']+=1
    else:raise AssertionError('Zero-alpha unavailable input did not raise: '+case+'/'+pair)
    out=analyze_bulk_scores(bulk,q,invalid_policy='complete_case',**call)
    selected=list(dict.fromkeys(sig[pair]['F']['genes']+sig[pair]['M']['genes']))
    mask=(bulk[selected]>0).all(axis=1)
    np.testing.assert_array_equal(out['valid_mask'].loc[bulk.index].to_numpy(),mask.to_numpy())
    assert out['coverage']['n_used']==int(mask.sum())
    assert set(out['exclusions'].sample_id)==set(bulk.index[~mask])
    assert out['scores'].loc[~mask].isna().all().all()
    rr,cr=compute_expected_on_mask(bulk,q,pair,alpha,method,mask)
    check(out['correlations']['raw']['r'],rr,'conditional_support');check(out['correlations']['corrected']['r'],cr,'conditional_support')
    counts['explicit_complete_case_calls']+=1
    continue
   out=analyze_bulk_scores(bulk,q,**call)
   assert out['coverage']['n_used']==118 and out['valid_mask'].all() and out['exclusions'].empty
   check(out['correlations']['raw']['r'],eraw.expected_r,'raw')
   check(out['correlations']['corrected']['r'],ecorr.expected_r,method)
   if pd.notna(eraw.archive_selector):check(out['correlations']['raw']['r'],eraw.archived_r,'archived_raw')
   if pd.notna(ecorr.archive_selector):check(out['correlations']['corrected']['r'],ecorr.archived_r,'archived_'+method)
   assert out['metadata']['leakage_status']=='disabled'
   counts['archive_comparable_calls']+=1

# Frozen Quah six-pair leakage fixture uses its own native compartment axis.
bulk=pd.read_csv(F/'Quah_paired_variable_variable/bulk_probability.csv',index_col=0,float_precision='round_trip').iloc[:3]
q=pd.DataFrame(np.tile([.2,.2,.6],(len(bulk),1)),index=bulk.index,columns=['Fibroblast','Macrophage','Epithelial'])
weighted=pd.read_csv(F/'expected_six_pair_weighted_indices.csv').set_index('pair_id')
lam=pd.read_csv(F/'expected_reference_gene_lambda.csv').set_index(['target','gene'])
for pair in sig:
 context={'case':'Quah_fixed_reference_q','pair_id':pair,'alpha':.01,'method':'known_loading'}
 out=analyze_bulk_scores(bulk,q,pair,bulk_scale='probability',alpha=.01,q_provenance='known',reference_law='product')
 assert out['metadata']['leakage_status'] not in ['not_estimated','disabled']
 gl=out['gene_leakage'];sl=out['score_leakage'];pl=out['pair_leakage']
 for _,row in gl.iterrows():
  context.update(sample_id=row.sample_id,signature=row.signature,gene=row.gene)
  oracle=lam.loc[(row.signature,row.gene)]
  check(row['lambda'],oracle.lambda_mean,'gene_lambda',2e-12)
  assert row.valid_states==oracle.n_valid_states and row.total_states==oracle.n_total_states
  counts['gene_lambda_rows_checked']+=1
 for _,row in sl.iterrows():
  context.update(sample_id=row.sample_id,signature=row.signature,gene='ALL_SCORE_GENES')
  w=weighted.loc[pair];side=row.signature
  check([row.weighted_sum,row.normalized_mean,row.beta],[w[f'weighted_lambda_{side}'],w[f'normalized_weighted_lambda_{side}'],w[f'beta_sum_{side}']],'weighted_index',2e-12)
 for _,row in pl.iterrows():check([row.weighted_pair_mean,row.normalized_pair_mean],[weighted.loc[pair,'weighted_lambda_pair'],weighted.loc[pair,'normalized_weighted_lambda_pair']],'weighted_index',2e-12)
 assert len(gl)==len(bulk)*sum(len(sig[pair][s]['genes']) for s in sig[pair])
 counts['six_pair_leakage_calls']+=1

# Independent permutation and unit-equivalence checks on a real safe P4 cohort.
bulk=pd.read_csv(F/'Quah_paired_variable_variable/bulk_probability.csv',index_col=0,float_precision='round_trip')
q=pd.read_csv(F/'Quah_paired_variable_variable/rna_q.csv',index_col=0,float_precision='round_trip')
for method in ['known_loading','joint_ols']:
 context={'case':'Quah_P4_unit_and_order','pair_id':'P4','alpha':0.,'method':method}
 base=analyze_bulk_scores(bulk,q,'P4',method=method,bulk_scale='probability',compute_leakage=False)
 changed=analyze_bulk_scores(bulk.iloc[::-1,::-1]*1e6,q.iloc[::-1],'P4',method=method,bulk_scale='cpm',compute_leakage=False)
 for key in ['raw','corrected']:check(base['correlations'][key]['r'],changed['correlations'][key]['r'],'unit_or_order_invariance')

implementation_hashes_after={str(p.relative_to(args.software_root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (args.software_root/'bulk_score_tools').rglob('*.py')}
assert implementation_hashes_before==implementation_hashes_after,'Implementation changed while validation was running; repeat on final code.'
result={'status':'PASS','validated_at_utc':datetime.now(timezone.utc).isoformat(),'counts':counts,'max_absolute_gaps':maxgap,'implementation_hashes':implementation_hashes_after,
 'scope':'Representative archived natural-profile simulation parity and frozen Quah reference diagnostics; not a new method-performance analysis.',
 'native_background_labels':'Quah Epithelial and Simkin Malignant are preserved. Correction parity disables reference leakage; Quah reference parity uses its own exact named axes.',
 'oracle_semantics':'known_q_OLS from small-alpha archive is only an oracle for joint_ols. Known-loading uses source logq_offset atalpha0 and independently computed fixed-weight subtraction atalpha.01.',
 'pure_target_boundary':'Some bulk scores are finite where the source pure log0 target was undefined. API numerical correlation availability does not assert pure-target identifiability.',
 'fixture_receipt_sha256':hashlib.sha256((F/'fixture_receipt.json').read_bytes()).hexdigest(),
 'API_SCHEMA_sha256':hashlib.sha256((args.software_root/'API_SCHEMA.md').read_bytes()).hexdigest()}
args.out.parent.mkdir(parents=True,exist_ok=True)
diff_path=args.out.with_name(args.out.stem+'_differences.csv');pd.DataFrame(diff_records).to_csv(diff_path,index=False,float_format='%.17g')
result['differences_csv']=diff_path.name;result['n_difference_checks']=len(diff_records)
args.out.write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result,indent=2))
