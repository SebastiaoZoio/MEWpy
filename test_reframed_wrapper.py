import unittest
import warnings
warnings.filterwarnings("ignore", message="scipy._lib.messagestream.MessageStream size changed")

from reframed import load_cbmodel

from mewpy.io import Reader, Engines, read_model
from mewpy.germ.analysis import FBA, fva, SRFBA, RFBA, single_gene_deletion, single_reaction_deletion

e_coli_gem_file = "/home/sebastiao_zoio/Documents/test_mewpy/mewpy/examples/models/germ/e_coli_core.xml"
e_coli_trn_model = "/home/sebastiao_zoio/Documents/test_mewpy/mewpy/examples/models/germ/e_coli_core_trn.csv"


def get_regulatory_reader(reg_file):
    return Reader(Engines.BooleanRegulatoryCSV,
                        reg_file,
                        sep=',',
                        id_col=0,
                        rule_col=2,
                        aliases_cols=[1],
                        header=0)


def get_orginal_gem_model(met_file, reg_file=None):
    gem_reader = Reader(Engines.MetabolicSBML, met_file)

    if reg_file:
        trn_reader = get_regulatory_reader(reg_file)
        return read_model(gem_reader, trn_reader)

    return read_model(gem_reader)



def get_new_gem_model(met_file, reg_file=None):
    reframed_model = load_cbmodel(met_file, use_infinity=False)
    reframed_reader = Reader(Engines.ReframedModelEngine, reframed_model)

    if reg_file:
        trn_reader = get_regulatory_reader(reg_file)
        return read_model(reframed_reader, trn_reader)

    return read_model(reframed_reader)


class ImportModelTest(unittest.TestCase):

    original_model = get_orginal_gem_model(e_coli_gem_file, e_coli_trn_model)
    new_model = get_new_gem_model(e_coli_gem_file, e_coli_trn_model)

    def setUp(self) -> None:
        pass

    def test_len_reactions(self):
        self.assertEqual(len(self.original_model.reactions),
                          len(self.new_model.reactions))
        
    def test_len_genes(self):
        self.assertEqual(len(self.original_model.genes),
                          len(self.new_model.genes))
        
    def test_len_metabolites(self):
        self.assertEqual(len(self.original_model.metabolites),
                          len(self.new_model.metabolites))
        
    def test_len_reg_metabolites(self):
        self.assertEqual(len(self.original_model.regulatory_metabolites),
                          len(self.new_model.regulatory_metabolites))
        
    def test_len_reg_reactions(self):
        self.assertEqual(len(self.original_model.regulatory_reactions),
                          len(self.new_model.regulatory_reactions))

        
    def test_equal_reactions(self):
        self.assertEqual(self.original_model.reactions.keys(), 
                         self.new_model.reactions.keys())

        for new_rxn in self.new_model.yield_reactions():
            original_rxn = self.original_model.reactions[new_rxn.id]

            self.assertEqual(original_rxn.id, new_rxn.id)
            self.assertEqual(original_rxn.types, new_rxn.types)
            self.assertEqual(original_rxn.equation, new_rxn.equation)
            self.assertEqual(original_rxn.bounds, new_rxn.bounds)
            self.assertEqual(original_rxn.reversibility, new_rxn.reversibility)
            self.assertEqual(original_rxn.metabolites.keys(), new_rxn.metabolites.keys())
            self.assertEqual(original_rxn.boundary, new_rxn.boundary)
            self.assertEqual(original_rxn.gpr.to_string(), new_rxn.gpr.to_string())
            self.assertEqual(original_rxn.genes.keys(), new_rxn.genes.keys())
            self.assertEqual(original_rxn.compartments, new_rxn.compartments)
            #self.assertEqual(original_rxn.charge_balance, new_rxn.charge_balance)
            self.assertEqual(original_rxn.mass_balance, new_rxn.mass_balance)



    def test_equal_metabolites(self):
        self.assertEqual(self.original_model.metabolites.keys(), 
                         self.new_model.metabolites.keys())

        for new_met in self.new_model.yield_metabolites():
            original_met = self.original_model.metabolites[new_met.id]
            
            self.assertEqual(original_met.id, new_met.id)
            self.assertEqual(original_met.types, new_met.types)
            self.assertEqual(original_met.compartment, new_met.compartment)
            self.assertEqual(original_met.formula, new_met.formula)
            self.assertEqual(original_met.molecular_weight, new_met.molecular_weight)
            #self.assertEqual(original_met.charge, new_met.charge)
            self.assertEqual(original_met.reactions.keys(), new_met.reactions.keys())



    def test_equal_genes(self):
        self.assertEqual(self.original_model.genes.keys(), 
                         self.new_model.genes.keys())

        for new_gene in self.new_model.yield_genes():
            original_gene = self.original_model.genes[new_gene.id]

            self.assertEqual(original_gene.id, new_gene.id)
            self.assertEqual(original_gene.types, new_gene.types)
            self.assertEqual(original_gene.coefficients, new_gene.coefficients)
            self.assertEqual(original_gene.is_active, new_gene.is_active)
            self.assertEqual(original_gene.reactions.keys(), new_gene.reactions.keys())

    
    def test_equal_reg_metabolites(self):
        self.assertEqual(self.original_model.regulatory_metabolites.keys(), 
                         self.new_model.regulatory_metabolites.keys())
        
        for key in self.new_model.regulatory_metabolites:
            original_met = self.original_model.regulatory_metabolites[key]
            new_met = self.new_model.regulatory_metabolites[key]

            self.assertEqual(original_met.id, new_met.id)
            self.assertEqual(original_met.interactions.keys(), new_met.interactions.keys())
            self.assertEqual(original_met.targets.keys(), new_met.targets.keys())


    def test_equal_reg_reactions(self):
        self.assertEqual(self.original_model.regulatory_reactions.keys(), 
                         self.new_model.regulatory_reactions.keys())
        
        for key in self.new_model.regulatory_reactions:
            original_rxn = self.original_model.regulatory_reactions[key]
            new_rxn = self.new_model.regulatory_reactions[key]

            self.assertEqual(original_rxn.id, new_rxn.id)
            self.assertEqual(original_rxn.coefficients, new_rxn.coefficients)
            self.assertEqual(original_rxn.is_active, new_rxn.is_active)
            self.assertEqual(original_rxn.interactions.keys(), new_rxn.interactions.keys())
            self.assertEqual(original_rxn.targets.keys(), new_rxn.targets.keys())
            self.assertEqual(original_rxn.environmental_stimulus, new_rxn.environmental_stimulus)

    

class MetabolicSimulationsTest(unittest.TestCase):
    
    original_model = get_orginal_gem_model(e_coli_gem_file, e_coli_trn_model)
    new_model = get_new_gem_model(e_coli_gem_file, e_coli_trn_model)

    def setUp(self) -> None:
        pass

    def test_fba(self):
        original_fba_problem = FBA(self.original_model).build()
        new_fba_problem = FBA(self.new_model).build()

        self.assertEqual(len(new_fba_problem.constraints), 0)
        self.assertEqual(len(new_fba_problem.variables), 0)

        original_fba_solution = original_fba_problem.optimize()
        new_fba_solution = new_fba_problem.optimize()

        self.assertAlmostEqual(original_fba_solution.objective_value, new_fba_solution.objective_value)

    
    def test_fva(self):
        original_fva = fva(self.original_model)
        new_fva = fva(self.new_model)

        for (new_index, new_row), (original_index, original_row) in zip(new_fva.iterrows(), original_fva.iterrows()):
            self.assertAlmostEqual(new_row["minimum"], original_row["minimum"])
            self.assertAlmostEqual(new_row["maximum"], original_row["maximum"])

    
    def test_srfba(self):
        original_fba_problem = SRFBA(self.original_model).build()
        new_fba_problem = SRFBA(self.new_model).build()

        self.assertEqual(len(original_fba_problem.variables), len(new_fba_problem.variables))
        self.assertEqual(len(original_fba_problem.constraints), len(new_fba_problem.constraints))

        original_fba_solution = original_fba_problem.optimize()
        new_fba_solution = new_fba_problem.optimize()
    
        self.assertAlmostEqual(original_fba_solution.objective_value, new_fba_solution.objective_value)


    def test_rfba(self):
        original_fba_problem = RFBA(self.original_model).build()
        new_fba_problem = RFBA(self.new_model).build()

        self.assertEqual(len(original_fba_problem.variables), len(new_fba_problem.variables))
        self.assertEqual(len(original_fba_problem.constraints), len(new_fba_problem.constraints))

        original_fba_solution = original_fba_problem.optimize()
        new_fba_solution = new_fba_problem.optimize()

        self.assertAlmostEqual(original_fba_solution.objective_value, new_fba_solution.objective_value)

        original_fba_dynamic_solution = original_fba_problem.optimize(dynamic=True)
        new_fba_dynamic_solution = new_fba_problem.optimize(dynamic=True)

        for key in original_fba_dynamic_solution.solutions:
            original_objective_value = original_fba_dynamic_solution.solutions[key].objective_value
            new_objective_value = new_fba_dynamic_solution.solutions[key].objective_value

        self.assertAlmostEqual(original_objective_value, new_objective_value)


class MetabolicAnalysisTest(unittest.TestCase):
    
    original_model = get_orginal_gem_model(e_coli_gem_file, e_coli_trn_model)
    new_model = get_new_gem_model(e_coli_gem_file, e_coli_trn_model)

    def setUp(self) -> None:
        pass


    def test_single_gene_deletion(self):
        new_sgd_df = single_gene_deletion(self.new_model)
        new_sgd_df['ids'] = new_sgd_df['ids'].str.join('')
        new_sgd_df['growth'] = new_sgd_df['growth'].fillna(0.0)
        new_sgd = new_sgd_df.set_index('ids').to_dict(orient='index')

        original_sgd_df = single_gene_deletion(self.original_model)
        original_sgd_df['status'] = original_sgd_df['status'].str.lower()
        original_sgd = original_sgd_df.to_dict(orient='index')  

        for gene, values in new_sgd.items():     
            self.assertAlmostEqual(values['growth'], original_sgd[gene]['growth'])
            self.assertEqual(values['status'], original_sgd[gene]['status'])


    def test_single_reaction_deletion(self):
        new_srd_df = single_reaction_deletion(self.new_model)
        new_srd_df['ids'] = new_srd_df['ids'].str.join('')
        new_srd_df['growth'] = new_srd_df['growth'].fillna(0.0)
        new_srd = new_srd_df.set_index('ids').to_dict(orient='index')

        original_srd_df = single_reaction_deletion(self.original_model)
        original_srd_df['status'] = original_srd_df['status'].str.lower()
        original_srd = original_srd_df.to_dict(orient='index')  

        for gene, values in new_srd.items():     
            self.assertAlmostEqual(values['growth'], original_srd[gene]['growth'])
            self.assertEqual(values['status'], original_srd[gene]['status'])

if __name__== '__main__':

    unittest.main()
