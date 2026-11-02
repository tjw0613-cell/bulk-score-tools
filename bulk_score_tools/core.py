"""No fitted scoring weights, hidden sampling, or selected-gene renormalization."""
from __future__ import annotations

from collections.abc import Mapping
from copy import deepcopy
from hashlib import sha256
from importlib.resources import files
import json
import math

import numpy as np
import pandas as pd


GENE_COLUMNS = ['sample_id','signature','target','gene','lambda','valid_states','total_states']
SCORE_COLUMNS = ['sample_id','signature','target','weighted_sum','normalized_mean',
                 'beta','defined_gene_count','total_gene_count','status']
PAIR_COLUMNS = ['sample_id','weighted_pair_mean','normalized_pair_mean','status']


def _labels(labels, name):
    labels = list(labels)
    ix = pd.Index(labels)
    if not labels or not ix.is_unique or ix.hasnans:
        raise ValueError(f'{name} must be nonempty, unique and nonmissing')
    if any(not isinstance(x, str) or not x.strip() for x in labels):
        raise ValueError(f'{name} must contain nonempty exact string labels')
    return labels


def _frame(frame, name):
    if not isinstance(frame,pd.DataFrame):
        raise ValueError(f'{name} must be a pandas DataFrame')
    if not frame.index.is_unique or frame.index.hasnans:
        raise ValueError(f'{name} sample IDs must be unique and nonmissing')
    _labels(frame.columns,f'{name} columns')
    try:
        x=frame.to_numpy(dtype=np.float64)
    except (ValueError,TypeError) as exc:
        raise ValueError(f'{name} must be numeric') from exc
    if not np.isfinite(x).all() or (x<0).any():
        raise ValueError(f'{name} must be finite and nonnegative')
    return x


def load_paper_signatures(preset):
    """Load pair-specific coefficients including the frozen Bill mean offset."""
    data=json.loads(files('bulk_score_tools').joinpath('data/paper_presets.json').read_text(encoding='utf-8'))
    if preset not in data['presets']:
        raise ValueError('Unknown paper preset; expected P1-P6')
    return deepcopy(data['presets'][preset]['signatures'])


def load_paper_reference():
    """Return the explicitly bundled 7-donor Quah full-human reference."""
    path=files('bulk_score_tools').joinpath('data/quah_reference.npz')
    with path.open('rb') as handle, np.load(handle,allow_pickle=False) as z:
        ref={key:z[key].copy() for key in z.files}
    for key in ['genes','donors','compartments']:
        ref[key]=ref[key].tolist()
    prov=json.loads(files('bulk_score_tools').joinpath('data/resource_provenance.json').read_text(encoding='utf-8'))
    ref.update(normalization='full_human_probability',version=prov['reference_version'],provenance=prov)
    return ref


def _resolve_signatures(signatures):
    preset=signatures if isinstance(signatures,str) else None
    if preset is not None:
        signatures=load_paper_signatures(preset)
    if not isinstance(signatures,Mapping) or len(signatures)!=2:
        raise ValueError('signatures must contain exactly two named definitions')
    _labels(signatures.keys(),'signature names')
    out={}
    for name,s in signatures.items():
        if not isinstance(s,Mapping) or 'target' not in s or 'genes' not in s:
            raise ValueError(f'{name} requires target and genes')
        target=_labels([s['target']],f'{name} target')[0]
        genes=_labels(s['genes'],f'{name} genes')
        w=s.get('weights',[1/len(genes)]*len(genes))
        if isinstance(w,Mapping):
            if set(w)!=set(genes):
                raise ValueError(f'{name} weight keys must match genes exactly')
            w=[w[g] for g in genes]
        try:
            w=np.asarray(w,dtype=float)
            offset=float(s.get('offset',0))
        except (ValueError,TypeError) as exc:
            raise ValueError(f'{name} has invalid coefficients') from exc
        if w.shape!=(len(genes),) or not np.isfinite(w).all() or (w<=0).any() or not np.isfinite(offset):
            raise ValueError(f'{name} weights must be finite and strictly positive; offset must be finite')
        with np.errstate(over='ignore',invalid='ignore'):
            beta=float(w.sum())
        if not math.isfinite(beta):
            raise ValueError(f'{name} coefficient sum beta must remain finite')
        out[name]={'target':target,'genes':genes,'weights':w.tolist(),'offset':offset,'beta':beta}
    return out,preset


def _pearson(a,b,n,reason=None,original_scales=None):
    if reason is None and n<3:
        reason='fewer_than_three_samples'
    if reason is None:
        # Rescale before centering: finite extreme coefficients must not overflow
        # means or Euclidean norms. The raw scale detects projection cancellation.
        scales=[float(np.max(np.abs(a))),float(np.max(np.abs(b)))]
        if original_scales is not None:
            scales=[max(float(original_scales[j]),scales[j]) for j in range(2)]
        aa=a/max(1.,scales[0]);bb=b/max(1.,scales[1])
        ac=aa-aa.mean();bc=bb-bb.mean()
        tol=64*np.finfo(float).eps
        if np.linalg.norm(ac)<=tol*max(1.,np.linalg.norm(aa)) or np.linalg.norm(bc)<=tol*max(1.,np.linalg.norm(bb)):
            reason='zero_numerical_variance'
    if reason is not None:
        return {'r':None,'status':'undefined','reason':reason,'n':int(n)}
    r=float(np.dot(ac,bc)/(np.linalg.norm(ac)*np.linalg.norm(bc)))
    return {'r':float(np.clip(r,-1,1)),'status':'ok','reason':None,'n':int(n)}


def _prepare_reference(reference,q_columns,genes,law,max_states):
    if not isinstance(reference,Mapping):
        raise ValueError('reference must be a versioned dictionary')
    required=['profiles','genes','donors','compartments','normalization','version','provenance']
    if any(k not in reference for k in required):
        raise ValueError('reference is missing required axes, normalization, version or provenance')
    if reference['normalization']!='full_human_probability':
        raise ValueError('reference must use the complete full-human probability denominator')
    if not isinstance(reference['version'],str) or not reference['version'].strip() or not isinstance(reference['provenance'],Mapping) or not reference['provenance']:
        raise ValueError('reference requires nonempty version and provenance')
    rg=_labels(reference['genes'],'reference genes')
    rc=_labels(reference['compartments'],'reference compartments')
    rd=_labels(reference['donors'],'reference donors')
    if set(rc)!=set(q_columns):
        raise ValueError('q and reference must have exactly the same compartment labels')
    if not set(genes)<=set(rg):
        raise ValueError('reference is missing selected genes; zero imputation is forbidden')
    x=np.asarray(reference['profiles'],dtype=float)
    if x.shape!=(len(rd),len(rc),len(rg)) or not np.isfinite(x).all() or (x<0).any() or (x.sum(axis=2)>1+1e-8).any():
        raise ValueError('invalid reference probability shape, values or full-human denominator')
    if law not in ('paired','product'):
        raise ValueError('reference_law must be paired or product')
    if isinstance(max_states,bool) or not isinstance(max_states,(int,np.integer)) or max_states<1:
        raise ValueError('max_reference_states must be a positive integer')
    ns=len(rd) if law=='paired' else len(rd)**len(rc)
    if ns>max_states:
        raise ValueError(f'exact reference enumeration requires {ns} states, exceeding {max_states}')
    if law=='paired':
        ix=np.repeat(np.arange(len(rd))[:,None],len(rc),axis=1)
    else:
        ix=np.indices((len(rd),)*len(rc)).reshape(len(rc),-1).T
    gi=[rg.index(g) for g in genes]
    component_profiles=x[ix,np.arange(len(rc))[None,:]][:,:,gi]
    h=sha256()
    h.update(np.ascontiguousarray(x,dtype='<f8').tobytes())
    h.update(json.dumps({'genes':rg,'donors':rd,'compartments':rc,'version':reference['version']},sort_keys=True).encode())
    meta={'version':reference['version'],'hash_sha256':h.hexdigest(),'law':law,
          'donor_count':len(rd),'compartments':rc,'state_count':ns,'provenance':deepcopy(dict(reference['provenance']))}
    return component_profiles,rc,meta


def _leakage(q,signatures,reference,law,max_states):
    genes=list(dict.fromkeys(g for s in signatures.values() for g in s['genes']))
    profiles,comps,meta=_prepare_reference(reference,list(q.columns),genes,law,max_states)
    q=q.loc[:,comps]
    gene_rows=[];score_rows=[];pair_rows=[]
    ns=profiles.shape[0]
    for sid,row in q.iterrows():
        contribution=profiles*row.to_numpy()[None,:,None]
        total=contribution.sum(axis=1)
        valid=total>0
        counts=valid.sum(axis=0)
        row_indices=[]
        for name,s in signatures.items():
            j=comps.index(s['target']);gj=[genes.index(g) for g in s['genes']]
            non=contribution[:,[k for k in range(len(comps)) if k!=j],:].sum(axis=1)
            ratios=np.divide(non,total,out=np.full_like(non,np.nan),where=valid)
            lam=np.divide(np.nansum(ratios,axis=0),counts,out=np.full(len(genes),np.nan),where=counts>0)
            vals=lam[gj];w=np.asarray(s['weights']);nvalid=int(np.isfinite(vals).sum())
            for g,gix in zip(s['genes'],gj):
                gene_rows.append([sid,name,s['target'],g,float(lam[gix]),int(counts[gix]),ns])
            good=nvalid==len(gj)
            weighted=float(np.dot(w,vals)) if good else float('nan')
            normalized=weighted/s['beta'] if good else float('nan')
            status='ok' if good else 'undefined_gene_leakage'
            score_rows.append([sid,name,s['target'],weighted,normalized,s['beta'],nvalid,len(gj),status])
            row_indices.append((weighted,normalized))
        ok=bool(np.isfinite(row_indices).all())
        pair_rows.append([sid,float(row_indices[0][0]/2+row_indices[1][0]/2) if ok else float('nan'),
                          float(row_indices[0][1]/2+row_indices[1][1]/2) if ok else float('nan'),
                          'ok' if ok else 'undefined_score_leakage'])
    return (pd.DataFrame(gene_rows,columns=GENE_COLUMNS),pd.DataFrame(score_rows,columns=SCORE_COLUMNS),
            pd.DataFrame(pair_rows,columns=PAIR_COLUMNS),meta)


def analyze_bulk_scores(bulk,q,signatures,*,reference=None,method='known_loading',bulk_scale='cpm',
                        alpha=0.0,full_human_totals=None,q_provenance='unspecified',q_scale='rna',
                        reference_law='product',invalid_policy='raise',compute_leakage=True,
                        max_reference_states=100000):
    """Score two fixed signatures and apply an explicitly chosen log-RNA correction.

    See API_SCHEMA.md for units, assumptions and the returned tables. Positive
    alpha defines a regularized score; it does not make zero-alpha recovery exact.
    """
    x=_frame(bulk,'bulk');_frame(q,'q')
    if method not in ('known_loading','joint_ols'):
        raise ValueError('method must be known_loading or joint_ols')
    if bulk_scale not in ('cpm','tpm','probability','counts'):
        raise ValueError('bulk_scale must be a linear cpm, tpm, probability or counts scale')
    if q_scale!='rna':
        raise ValueError('q_scale must be rna; cell fractions require an upstream RNA-content conversion')
    if q_provenance not in ('known','estimated','unspecified'):
        raise ValueError('q_provenance must be known, estimated or unspecified')
    if invalid_policy not in ('raise','complete_case'):
        raise ValueError('invalid_policy must be raise or complete_case')
    if not isinstance(compute_leakage,(bool,np.bool_)):
        raise ValueError('compute_leakage must be boolean')
    if reference_law not in ('paired','product'):
        raise ValueError('reference_law must be paired or product')
    try:
        alpha=float(alpha)
    except (ValueError,TypeError) as exc:
        raise ValueError('alpha must be finite and nonnegative') from exc
    if not math.isfinite(alpha) or alpha<0:
        raise ValueError('alpha must be finite and nonnegative')
    s,preset=_resolve_signatures(signatures);names=list(s)
    if len(q)!=len(bulk) or not bulk.index.isin(q.index).all():
        raise ValueError('q and bulk must have the same sample IDs')
    q=q.loc[bulk.index].copy()
    qv=q.to_numpy(dtype=float)
    if not np.allclose(qv.sum(axis=1),1,atol=1e-8,rtol=0):
        raise ValueError('complete RNA fractions must sum to 1; no silent renormalization')
    targets=[s[name]['target'] for name in names]
    if not set(targets)<=set(q.columns):
        raise ValueError('q is missing a target compartment')
    selected=list(dict.fromkeys(g for v in s.values() for g in v['genes']))
    if not set(selected)<=set(bulk.columns):
        raise ValueError('bulk is missing selected genes; zero imputation is forbidden')
    normalization='caller_supplied_complete_denominator'
    if bulk_scale=='counts':
        if not isinstance(full_human_totals,pd.Series) or not full_human_totals.index.is_unique or full_human_totals.index.hasnans:
            raise ValueError('counts require a labelled full_human_totals Series')
        if len(full_human_totals)!=len(bulk) or not bulk.index.isin(full_human_totals.index).all():
            raise ValueError('full_human_totals must match the sample IDs')
        total=full_human_totals.loc[bulk.index].to_numpy(dtype=float)
        if not np.isfinite(total).all() or (total<=0).any() or (x.sum(axis=1)>total*(1+1e-8)).any():
            raise ValueError('invalid full-human counts totals')
        x=x/total[:,None]*1e6
        normalization='supplied_full_human_totals'
    elif full_human_totals is not None:
        raise ValueError('full_human_totals is only used with counts')
    elif bulk_scale=='probability':
        x=x*1e6
    if (x.sum(axis=1)>1e6*(1+1e-8)).any():
        raise ValueError('selected linear expression exceeds the complete per-million total')
    reasons=[[] for _ in range(len(bulk))]
    if alpha==0:
        for i,g in enumerate(bulk.columns):
            if g in selected:
                for r in np.flatnonzero(x[:,i]==0):reasons[r].append(f'zero_expression:{g}')
    for t in dict.fromkeys(targets):
        for r in np.flatnonzero(q[t].to_numpy()<=0):reasons[r].append(f'nonpositive_target_q:{t}')
    valid=np.array([not r for r in reasons],dtype=bool)
    excluded=pd.DataFrame([[sid,';'.join(r)] for sid,r in zip(bulk.index,reasons) if r],columns=['sample_id','reasons'])
    if not valid.all() and invalid_policy=='raise':
        raise ValueError(f'unscorable samples: {excluded.to_dict(orient="records")}')
    scores=pd.DataFrame(np.nan,index=bulk.index,columns=[f'{name}__{kind}' for name in names for kind in ['raw','corrected']])
    usable=x[valid];q_used=q.loc[valid];n=int(valid.sum())
    raw=[]
    for name in names:
        v=s[name];idx=[bulk.columns.get_loc(g) for g in v['genes']]
        with np.errstate(over='ignore',invalid='ignore'):
            a=np.log2(usable[:,idx]+alpha)@np.asarray(v['weights'])+v['offset']
        if not np.isfinite(a).all():
            raise ValueError('nonfinite derived raw score; check coefficient and expression scales')
        raw.append(a)
        scores.loc[valid,f'{name}__raw']=a
    raw=np.column_stack(raw)
    logq=np.log2(q_used.loc[:,targets].to_numpy(dtype=float))
    ols=None;corr_reason=None
    if method=='known_loading':
        with np.errstate(over='ignore',invalid='ignore'):
            adjusted=raw-logq*np.array([s[name]['beta'] for name in names])[None,:]
    else:
        design=np.column_stack([np.ones(n),logq])
        if n:
            coef,_,rank,_=np.linalg.lstsq(design,raw,rcond=None)
            with np.errstate(over='ignore',invalid='ignore'):
                adjusted=raw-design@coef
        else:
            rank=0;adjusted=raw.copy()
        df=n-int(rank)
        with np.errstate(over='ignore',invalid='ignore'):
            residual_var={name:float(np.var(adjusted[:,j],ddof=1)) if n>1 else None for j,name in enumerate(names)}
        if any(v is not None and not math.isfinite(v) for v in residual_var.values()):
            raise ValueError('nonfinite derived residual variance; check numerical scales')
        ols={'design_columns':['intercept',*[f'log2_q:{t}' for t in targets]],'rank':int(rank),
             'residual_df':df,'rank_deficient':bool(rank<3),
             'residual_variance':residual_var}
        if df<=0:corr_reason='nonpositive_residual_degrees_of_freedom'
    if not np.isfinite(adjusted).all():
        raise ValueError('nonfinite derived corrected score; check coefficient and composition scales')
    for j,name in enumerate(names):scores.loc[valid,f'{name}__corrected']=adjusted[:,j]
    correlations={'raw':_pearson(raw[:,0],raw[:,1],n),
                  'corrected':_pearson(adjusted[:,0],adjusted[:,1],n,corr_reason,
                                       original_scales=np.max(np.abs(raw),axis=0) if n else np.zeros(2))}
    gene=pd.DataFrame(columns=GENE_COLUMNS);score=pd.DataFrame(columns=SCORE_COLUMNS);pair=pd.DataFrame(columns=PAIR_COLUMNS)
    refmeta=None;leakstatus='not_estimated'
    if not compute_leakage:
        leakstatus='disabled'
    else:
        if reference is None and preset is not None:reference=load_paper_reference()
        if reference is not None:
            gene,score,pair,refmeta=_leakage(q_used,s,reference,reference_law,max_reference_states)
            leakstatus='reference_conditioned'
    meta={'method':method,'alpha':alpha,'bulk_scale':bulk_scale,'normalization_source':normalization,
          'q_scale':q_scale,'q_provenance':q_provenance,'signature_order':names,'signatures':s,
          'shared_genes':sorted(set(s[names[0]]['genes'])&set(s[names[1]]['genes'])),
          'preset':preset,'leakage_status':leakstatus,'reference':refmeta,'ols':ols,
          'q_max_sum_deviation':float(np.max(np.abs(qv.sum(axis=1)-1))) if len(qv) else 0.,
          'interpretation':[
              'Numeric q and observed bulk do not establish the conditions for exact pure-score recovery.',
              'Known-loading correction preserves supplied fixed coefficients; joint OLS also removes state variation associated with composition.',
              'Positive alpha defines a regularized score target; no zero-alpha guarantee is asserted.',
              'Reference leakage is descriptive and reference-conditioned, not patient-specific truth, risk classification or calibrated recovery error.']}
    return {'schema_version':'1.0','correlations':correlations,'scores':scores,
            'valid_mask':pd.Series(valid,index=bulk.index,name='valid'), 'exclusions':excluded,
            'coverage':{'n_input':len(bulk),'n_used':n,'n_excluded':len(bulk)-n,
                        'selected_gene_count':len(selected),'signature_gene_counts':{k:len(v['genes']) for k,v in s.items()}},
            'gene_leakage':gene,'score_leakage':score,'pair_leakage':pair,'metadata':meta}
