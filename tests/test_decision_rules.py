"""Synthetic engineering fixtures; never evidence of biomedical performance."""
import json
from pathlib import Path
import tempfile
import unittest
import uuid
import shutil
from unittest.mock import patch
import decision_rules as rules
from constrained_learning import train
from scientist_workflow import review

def rule(kind='gene_role', value='injury_marker', subject='GENE1', **extra):
    return dict(kind=kind,value=value,subject=subject,disease='TEST',species='human',tissue='blood',stage='early',
                source_id='SYNTHETIC_TEST',quote='synthetic source quotation',source_text='synthetic source quotation',**extra)

CONTEXT=dict(disease='TEST',species='human',tissue='blood',stage='early')

class Decisions(unittest.TestCase):
    def priority(self):
        return {'filters':{},'eligible_genes':[{'gene':'GENE1','tier':'medium','sources':[]}],'high_priority_genes':[], 'medium_priority_genes':[]}

    def test_context_unknown_does_not_exclude(self):
        out=rules.apply_gene_rules(self.priority(),[rule()],{'disease':'TEST'})
        self.assertEqual(len(out['eligible_genes']),1)
        self.assertEqual(out['eligible_genes'][0]['scientific_decision']['role'],'unknown')

    def test_held_cannot_reenter_drug_or_herb_lookup(self):
        from gene_drug_discovery import collect_gene_evidence
        from target_based_screen import reverse_lookup_genes
        out=rules.apply_gene_rules(self.priority(),[rule()],CONTEXT)
        self.assertEqual(out['drug_lookup_gene_count'],0)
        self.assertEqual(len(out['scientifically_held_genes']),1)
        report={'gene_evidence_prioritization':out}
        self.assertEqual(collect_gene_evidence(report),{})
        self.assertEqual(reverse_lookup_genes(report),[])

    def test_conflicts_and_direction_mismatch(self):
        out=rules.apply_gene_rules(self.priority(),[rule(),rule(value='driver')],CONTEXT)
        self.assertEqual(len(out['eligible_genes']),0)
        out=rules.apply_gene_rules(self.priority(),[rule('target_direction','inhibit')],CONTEXT)
        drug=rules.apply_drug_rules({'candidates':[{'gene':'GENE1','drug':'DRUG1'}]},out,
             [rule('drug_direction','activate','DRUG1',target_gene='GENE1')],CONTEXT)
        self.assertEqual(drug['candidates'][0]['scientific_decision']['state'],'held')

    def test_no_rules_preserves_candidates_but_exposes_unknowns(self):
        priority=rules.apply_gene_rules(self.priority(),[],CONTEXT)
        out=rules.apply_drug_rules({'candidates':[{'gene':'GENE1','drug':'DRUG1'}]},priority,[],CONTEXT)
        self.assertEqual(out['candidates'][0]['scientific_decision']['state'],'requires_evidence')

    def test_exposure_rejects_wrong_scale_and_nonfinite(self):
        r=rule('exposure','insufficient','DRUG1',comparison_basis='total_plasma_vs_kd')
        with self.assertRaises(ValueError):rules.validate(r)
        r.update(comparison_basis='unbound_site_vs_unbound_effective',dose='1',time='1h',assay='matched',site_unbound_nM=float('nan'),effective_unbound_nM=2)
        with self.assertRaises(ValueError):rules.validate(r)

    def test_versions_are_immutable_and_owner_isolated(self):
        folder=Path(__file__).resolve().parent / ('rules-test-' + uuid.uuid4().hex)
        folder.mkdir()
        self.addCleanup(shutil.rmtree, folder)
        with patch.object(rules,'DB',folder/'rules.sqlite'):
            old=rules.save_revision('a',[rule()],reason='first')
            new=rules.save_revision('a',[],parent=old['id'],reason='correction')
            self.assertEqual(len(rules.get_revision('a',old['id'])['rules']),1)
            self.assertEqual(new['parent'],old['id'])
            with self.assertRaises(ValueError):rules.get_revision('b',old['id'])

    def test_three_roles_cover_all_records(self):
        calls=[]
        def model(cfg,messages,**kw):
            entries=json.loads(messages[0]['content'])['records'];calls.append(kw['system'])
            return json.dumps({'reviews':[{'ids':[r['id'] for r in entries], 'hypothesis':'synthetic',
                     'counterevidence_or_unknown':'unknown','next_test':'independent test'}]})
        out=review('gene','TEST',[{'gene':str(i)} for i in range(23)],{},model)
        self.assertEqual(out['processed_count'],23)
        self.assertEqual(len(calls),9)
        self.assertEqual(len(out['input_records']),23)
        self.assertEqual(len({i['id'] for i in out['input_records']}),23)

    def test_offline_training_and_leakage_rejection(self):
        rows=[{'label':i%2,'adjudicated':True,'source_group':f'g{i//8}', 'source_id':f's{i//8}',
               'features':{'direction_match':float(i%2),'exposure_match':None}} for i in range(80)]
        out=train({'task':'drug','rows':rows})
        self.assertTrue(out['fit_performed'])
        self.assertIn('logistic_baseline_validation',out)
        self.assertTrue(all(v>=0 for row in out['weights']['w1'][:6] for v in row))
        rows[1]['source_group']='different'
        with self.assertRaises(ValueError):train({'task':'drug','rows':rows})

if __name__=='__main__':unittest.main()
