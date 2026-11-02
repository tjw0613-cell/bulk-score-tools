"""Run the bundled deterministic synthetic P4 example, with no external files."""
from pathlib import Path
import argparse
import json
import pandas as pd
from bulk_score_tools import analyze_bulk_scores


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path(__file__).resolve().parent/'example_output')
    parser.add_argument('--alpha',type=float,default=.01,help='Pseudocount in CPM units; demo default .01')
    parser.add_argument('--method',choices=['known_loading','joint_ols'],default='known_loading')
    args=parser.parse_args()
    root=Path(__file__).resolve().parent
    bulk=pd.read_csv(root/'examples/bulk_cpm.csv',index_col='sample_id')
    q=pd.read_csv(root/'examples/q_rna.csv',index_col='sample_id')
    result=analyze_bulk_scores(bulk,q,'P4',method=args.method,alpha=args.alpha,
                               bulk_scale='cpm',q_provenance='known')
    args.output.mkdir(parents=True,exist_ok=True)
    result['scores'].to_csv(args.output/'scores.csv',index_label='sample_id',float_format='%.17g')
    result['valid_mask'].to_csv(args.output/'valid_mask.csv',index_label='sample_id')
    for key in ['gene_leakage','score_leakage','pair_leakage','exclusions']:
        result[key].to_csv(args.output/f'{key}.csv',index=False,float_format='%.17g')
    summary={key:result[key] for key in ['schema_version','correlations','coverage','metadata']}
    (args.output/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({'example':'synthetic_reference_mixtures_not_clinical_validation',
                      'alpha':args.alpha,'method':args.method,'correlations':result['correlations'],
                      'output':str(args.output)},ensure_ascii=True,indent=2))


if __name__=='__main__':main()
