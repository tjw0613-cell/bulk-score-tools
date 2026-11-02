import copy
import unittest
import numpy as np
import pandas as pd
from numpy.testing import assert_allclose
from bulk_score_tools import analyze_bulk_scores, load_paper_reference, load_paper_signatures


def toy():
    ids = [f's{i}' for i in range(8)]
    q = pd.DataFrame({'F': [.12,.17,.22,.28,.31,.35,.39,.44],
                      'M': [.43,.35,.31,.27,.22,.18,.15,.12]}, index=ids)
    q['B'] = 1-q.sum(axis=1)
    pure = pd.DataFrame({'x': [2,4,3,9,5,8,11,7], 'y': [7,4,8,2,5,6,3,10]}, index=ids)
    bulk = pure.mul(q[['F','M']].set_axis(['x','y'],axis=1))
    sig = {'sf': {'target':'F','genes':['x'],'weights':[.7],'offset':-1},
           'sm': {'target':'M','genes':['y'],'weights':[1.3],'offset':2}}
    return bulk, q, sig, pure


def ref():
    # At q=(.5,.5), ratios are 1/2, 9/10, 1/10, 1/2: mean 1/2.
    # Ratio of mean expressions is also 1/2 globally, so use unequal x states below.
    return {'profiles': np.array([[[1e-6,0],[1e-6,2e-6]],
                                  [[9e-6,0],[3e-6,4e-6]]]),
            'genes':['x','y'], 'donors':['d1','d2'], 'compartments':['F','M'],
            'normalization':'full_human_probability','version':'gold-1',
            'provenance':{'source':'deterministic test fixture'}}


class ScoreTests(unittest.TestCase):
    def test_known_loading_pointwise_with_state_composition_dependence(self):
        b,q,s,p = toy()
        r=analyze_bulk_scores(b,q,s)
        assert_allclose(r['scores']['sf__corrected'], .7*np.log2(p.x)-1, atol=2e-14)
        assert_allclose(r['scores']['sm__corrected'], 1.3*np.log2(p.y)+2, atol=2e-14)
        self.assertEqual(r['metadata']['leakage_status'],'not_estimated')

    def test_units_full_denominator(self):
        b,q,s,_=toy(); totals=pd.Series(np.arange(8)*1234+100000,index=b.index)
        base=analyze_bulk_scores(b,q,s)
        prob=analyze_bulk_scores(b/1e6,q,s,bulk_scale='probability')
        count=analyze_bulk_scores(b.mul(totals,axis=0)/1e6,q,s,bulk_scale='counts',full_human_totals=totals)
        assert_allclose(base['scores'],prob['scores'],atol=2e-14)
        assert_allclose(base['scores'],count['scores'],atol=2e-14)
        with self.assertRaises(ValueError): analyze_bulk_scores(b,q,s,bulk_scale='counts')

    def test_permutation_and_weight_scaling(self):
        b,q,s,_=toy(); r=analyze_bulk_scores(b,q,s)
        z=analyze_bulk_scores(b.iloc[::-1,::-1],q.sample(frac=1,random_state=1),s)
        assert_allclose(z['scores'].loc[b.index],r['scores'],atol=2e-14)
        t=copy.deepcopy(s)
        for v in t.values(): v['weights']=[x*3 for x in v['weights']]; v['offset']*=3
        zz=analyze_bulk_scores(b,q,t)
        assert_allclose(zz['scores'],r['scores']*3,atol=2e-14)
        self.assertAlmostEqual(zz['correlations']['corrected']['r'],r['correlations']['corrected']['r'])

    def test_alpha_is_per_million(self):
        b,q,s,_=toy(); b.iloc[0,0]=0
        with self.assertRaises(ValueError): analyze_bulk_scores(b,q,s)
        r=analyze_bulk_scores(b,q,s,alpha=.01)
        self.assertAlmostEqual(r['scores'].iloc[0,0],.7*np.log2(.01)-1)
        a=analyze_bulk_scores(b/1e6,q,s,alpha=.01,bulk_scale='probability')
        assert_allclose(r['scores'],a['scores'],atol=2e-14)

    def test_common_mask_and_reasons(self):
        b,q,s,_=toy(); b.iloc[0,0]=0; q.loc['s1',['F','B']]=[0,.65]
        r=analyze_bulk_scores(b,q,s,invalid_policy='complete_case')
        self.assertEqual(r['coverage']['n_used'],6)
        self.assertEqual(r['exclusions'].sample_id.tolist(),['s0','s1'])
        self.assertTrue(r['scores'].iloc[:2].isna().all().all())

    def test_fixed_q_ols_and_degenerate_variance(self):
        b,q,s,_=toy(); q.loc[:]=[.2,.2,.6]
        r=analyze_bulk_scores(b,q,s,method='joint_ols')
        self.assertEqual(r['metadata']['ols']['rank'],1)
        self.assertAlmostEqual(r['correlations']['raw']['r'],r['correlations']['corrected']['r'])
        b.loc[:]=1
        r=analyze_bulk_scores(b,q,s,method='joint_ols')
        self.assertIsNone(r['correlations']['corrected']['r'])
        self.assertEqual(r['correlations']['corrected']['status'],'undefined')

    def test_ols_does_not_equal_known_loading(self):
        b,q,s,p=toy()
        k=analyze_bulk_scores(b,q,s)
        o=analyze_bulk_scores(b,q,s,method='joint_ols')
        self.assertGreater(abs(k['correlations']['corrected']['r']-o['correlations']['corrected']['r']),.01)

    def test_exact_ols_span_large_coefficients_is_undefined(self):
        b,q,s,_=toy()
        b['x']=4*q.F;b['y']=8*q.M
        s['sf']['weights']=[1e5];s['sm']['weights']=[2e5]
        s['sf']['offset']=1e7;s['sm']['offset']=-2e7
        r=analyze_bulk_scores(b,q,s,method='joint_ols')
        self.assertIsNone(r['correlations']['corrected']['r'])
        self.assertEqual(r['correlations']['corrected']['reason'],'zero_numerical_variance')

    def test_small_resolved_variance_remains_defined(self):
        b,q,s,p=toy();q.loc[:]=[.2,.2,.6]
        b['x']=np.exp2(np.arange(8)*1e-8)
        b['y']=np.exp2(np.array([2,4,1,5,3,7,6,0])*1e-8)
        r=analyze_bulk_scores(b,q,s,method='joint_ols')
        self.assertEqual(r['correlations']['corrected']['status'],'ok')

    def test_nonfinite_derived_coefficients_and_scores_fail(self):
        b,q,s,_=toy()
        t=copy.deepcopy(s);t['sf'].update(genes=['x','y'],weights=[1e308,1e308])
        with self.assertRaises(ValueError):analyze_bulk_scores(b,q,t)
        t=copy.deepcopy(s);t['sf']['weights']=[1e308]
        with self.assertRaises(ValueError):analyze_bulk_scores(b*100,q,t)
        b['x']=1
        t['sf']['offset']=1e308
        with self.assertRaises(ValueError):analyze_bulk_scores(b,q,t)

    def test_large_finite_scaling_keeps_correlation_defined(self):
        b,q,s,_=toy();r=analyze_bulk_scores(b,q,s)
        t=copy.deepcopy(s)
        for v in t.values():v['weights']=[w*1e200 for w in v['weights']];v['offset']*=1e200
        z=analyze_bulk_scores(b,q,t)
        self.assertEqual(z['correlations']['corrected']['status'],'ok')
        self.assertAlmostEqual(z['correlations']['corrected']['r'],r['correlations']['corrected']['r'])

    def test_bad_inputs_fail(self):
        b,q,s,_=toy()
        cases=[{'q_scale':'cell'}, {'alpha':-1}, {'bulk_scale':'log'}, {'q_provenance':'pretend_true'}, {'method':'unknown'}]
        for kw in cases:
            with self.subTest(kw=kw),self.assertRaises(ValueError): analyze_bulk_scores(b,q,s,**kw)
        for bad in [b.set_axis(['z']*8), b.rename(columns={'y':'x'}), b.drop(columns='x'), b.assign(x=-1), b.assign(x=np.inf)]:
            with self.subTest(bad='bulk'),self.assertRaises(ValueError): analyze_bulk_scores(bad,q,s)
        badq=q.copy(); badq.F*=2
        with self.assertRaises(ValueError): analyze_bulk_scores(b,badq,s)
        badq=q.copy(); badq.index=[*q.index[:-1],'other']
        with self.assertRaises(ValueError): analyze_bulk_scores(b,badq,s)

    def test_too_few_samples_is_undefined(self):
        b,q,s,_=toy();r=analyze_bulk_scores(b.iloc[:2],q.iloc[:2],s)
        self.assertEqual(r['correlations']['raw']['reason'],'fewer_than_three_samples')


class LeakageTests(unittest.TestCase):
    def call(self,reference=None,**kw):
        b=pd.DataFrame({'x':[1,2,3], 'y':[4,2,5]},index=['a','b','c'])
        q=pd.DataFrame({'F':[.5]*3,'M':[.5]*3},index=b.index)
        s={'f':{'target':'F','genes':['x']},'m':{'target':'M','genes':['y']}}
        return analyze_bulk_scores(b,q,s,reference=ref() if reference is None else reference,**kw)

    def test_ratio_before_mean_gold(self):
        r=self.call();g=r['gene_leakage'];x=g[(g.sample_id=='a')&(g.gene=='x')].iloc[0]
        want=np.mean([1/2,3/4,1/10,3/12])
        self.assertAlmostEqual(x['lambda'],want)
        self.assertNotAlmostEqual(x['lambda'],2/(5+2))
        self.assertEqual(x.valid_states,4)
        self.assertTrue(g['lambda'].between(0,1).all())

    def test_paired_law_gold(self):
        r=self.call(reference_law='paired');g=r['gene_leakage']
        self.assertAlmostEqual(g[(g.sample_id=='a')&(g.gene=='x')]['lambda'].item(),.375)

    def test_undefined_states_not_zero_filled(self):
        rr=ref();rr['profiles'][:,:,0]=0
        r=self.call(reference=rr); x=r['gene_leakage'].query("gene=='x'")
        self.assertTrue(x['lambda'].isna().all());self.assertTrue((x.valid_states==0).all())
        self.assertTrue(r['pair_leakage'].weighted_pair_mean.isna().all())

    def test_partial_defined_states_reported(self):
        rr=ref();rr['profiles'][0,:,0]=0
        r=self.call(reference=rr);x=r['gene_leakage'].query("gene=='x'").iloc[0]
        self.assertEqual(x.valid_states,3)
        self.assertAlmostEqual(x['lambda'],(1+0+3/12)/3)

    def test_reference_cost_and_missing_features_fail(self):
        with self.assertRaises(ValueError):self.call(max_reference_states=3)
        rr=ref();rr['genes']=['z','y']
        with self.assertRaises(ValueError):self.call(reference=rr)

    def test_shared_genes_reported_unchanged(self):
        b,q,s,_=toy();s['sm']['genes']=['x'];s['sm']['weights']=[1.3]
        r=analyze_bulk_scores(b,q,s)
        self.assertEqual(r['metadata']['shared_genes'],['x'])

    def test_presets_are_fixed_resources(self):
        rr=load_paper_reference();s=load_paper_signatures('P4')
        self.assertEqual(rr['profiles'].shape,(7,3,33))
        self.assertEqual(len(s['F']['genes']),8)
        self.assertTrue(np.all(np.array(s['F']['weights'])>0))
        b=pd.DataFrame(np.ones((3,33)),index=['a','b','c'],columns=rr['genes'])
        q=pd.DataFrame([[.2,.2,.6]]*3,index=b.index,columns=rr['compartments'])
        r=analyze_bulk_scores(b,q,'P4',alpha=.01)
        self.assertEqual(r['metadata']['reference']['state_count'],343)
        self.assertEqual(r['metadata']['preset'],'P4')


if __name__=='__main__':unittest.main(verbosity=2)
