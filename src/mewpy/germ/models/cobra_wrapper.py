from typing import Any, Sequence, TYPE_CHECKING, Generator, Union, List

from .metabolic_wrapper import MetabolicModelWrapper
from mewpy.util.utilities import generator
from mewpy.util.history import HistoryManager
from mewpy.germ.algebra import Expression
from mewpy.germ.variables.variable import Variable
from mewpy.germ.models.serialization import serialize
from mewpy.simulation.cobra import Simulation
from mewpy.io.engines.engines_utils import build_symbolic
from mewpy.solvers.solution import Solution

if TYPE_CHECKING:
    from mewpy.germ.variables import Gene, Metabolite, Reaction


class CobraModelWrapper(MetabolicModelWrapper, model_type='metabolic_wrapper', register=True, constructor=True, checker=True):
    """
    A germ metabolic wrapper model consists of a reference to a classic Genome-Scale Metabolic (GEM) model,
    represented as a COBRApy or Reframed object.

    GEM models are systems biology tools used to predict the phenotype of an organism or cellular community
    in range of environmental and genetic conditions.
    To perform phenotype prediction, a metabolic model can be attached to several simulation methods:
        - FBA
        - pFBA
        - FVA
        - ...


    The metabolic model, as with other models, provides a clean interface for manipulation with the add, remove and
    update methods. One can perform the following operations:
        - Add reactions, metabolites and genes
        - Remove reactions, metabolites and genes
        - Update the objective function
    """
    def __init__(self,
                 identifier: Any,
                 **kwargs):

        """
        A germ metabolic model consists of a classic Genome-Scale Metabolic (GEM) model,
        containing reactions, metabolites and genes.

        GEM models are systems biology tools used to predict the phenotype of an organism or cellular community
        in range of environmental and genetic conditions.
        To perform phenotype prediction, a metabolic model can be attached to several simulation methods:
            - FBA
            - pFBA
            - FVA
            - ...


        The metabolic model, as with other models, provides a clean interface for manipulation with the add, remove and
        update methods. One can perform the following operations:
            - Add reactions, metabolites and genes
            - Remove reactions, metabolites and genes
            - Update the objective function

        :param identifier: identifier, e.g. iMC1010
        :param metabolic_reference: a reference to an instance of a metabolic network object instance
        """


        super().__init__(identifier,
                         **kwargs)
        
        self._met_model = None
        self._genes = {}
        self._metabolites = {}
        self._reactions = {}

        self._parsed = 0


    def set_simulator(self):
        self.simulator = Simulation(self._met_model)

    # -----------------------------------------------------------------------------
    # Model type manager
    # -----------------------------------------------------------------------------
    @serialize('types', None)
    @property
    def types(self):
        """
        Returns the types of the model
        :return: a set with the types of the model
        """
        _types = {MetabolicModelWrapper.model_type}

        _types.update(super(MetabolicModelWrapper, self).types)

        return _types


    @property
    def parsed(self):
        try:
            is_parsed = self._parsed
        except:
            is_parsed = False
        
        return is_parsed
        
    @property
    def met_model(self):
        return self.simulator.model

    @property
    def metabolites(self):
        if not self._metabolites:
            self._metabolites = {cobra_met: self.cobra_met_to_mewpy_met(cobra_met) for cobra_met in self.simulator.metabolites}
        return self._metabolites
    
    @property
    def genes(self):
        if not self._genes:
            self._genes = {cobra_gene: self.cobra_gene_to_mewpy_gene(cobra_gene) for cobra_gene in self.simulator.genes}
        return self._genes

    @property
    def reactions(self):
        if not self._reactions:
            self._reactions = {cobra_rxn: self.cobra_rxn_to_mewpy_rxn(cobra_rxn) for cobra_rxn in self.simulator.reactions}
        return self._reactions

    

    @property
    def objective(self):
        return self.cobra_objective_to_objective()
    
    @property
    def contexts(self) -> List[HistoryManager]:
        return super(MetabolicModelWrapper, self).contexts
    

    def yield_genes(self) -> Generator['Gene', None, None]:
        """
        It yields the genes of the model.
        :return: a generator with the genes of the model
        """
        return generator(self.genes)

    def yield_gprs(self) -> Generator['Expression', None, None]:
        """
        It yields the GPRs of the model.
        :return: a generator with the GPRs of the model
        """
        return (value.gpr for value in self.reactions.values())

    def yield_metabolites(self) -> Generator['Metabolite', None, None]:
        """
        It yields the metabolites of the model.
        :return: a generator with the metabolites of the model
        """
        return generator(self.metabolites)

    def yield_reactions(self) -> Generator['Reaction', None, None]:
        """
        It yields the reactions of the model.
        :return: a generator with the reactions of the model
        """
        return generator(self.reactions)


    def add_reg_data(self, args, id):
        if self.is_regulatory():
            if id in self.regulators:
                args['types'].add('regulator')
                args['interactions'] = self.regulators[id].interactions
    

            if id in self.targets:
                args['types'].add('target')
                args['interaction'] = self.targets[id].interaction


    def add_reaction_data(self, args:dict, rxn_id):
        rxn_dict = self.simulator.get_reaction(rxn_id)

        symbolic, warning = build_symbolic(expression=rxn_dict['gpr'])
        genes = {gene.name: self.cobra_gene_to_mewpy_gene(gene.name) for gene in symbolic.atoms(symbols_only=True)}
        gpr = Expression(symbolic=symbolic, variables=genes)

        stoichiometry = {self.cobra_met_to_mewpy_met(met_id): coeff for met_id, coeff in rxn_dict['stoichiometry'].items()}

        args['types']=set(['reaction'])
        args['identifier']=rxn_id
        args['name']=rxn_dict['name']
        args['aliases']={rxn_id, rxn_dict['name']}
        args['model']=self
        args['bounds']=(rxn_dict['lb'], rxn_dict['ub'])
        args['stoichiometry']=stoichiometry
        args['gpr']=gpr

            
    def add_metabolite_data(self, args:dict, met_id):
        met_dict = self.simulator.get_metabolite(met_id)

        reactions = {rxn: self.mewpy_rxn_instance(rxn) for rxn in self.simulator.get_metabolite_reactions(met_id)}

        args['types'] = set(['metabolite'])
        args['identifier'] = met_id
        args['name'] = met_dict['name']
        args['aliases'] = {met_id, met_dict['name']}
        args['model'] = self
        args['formula'] = met_dict['formula']
        args['compartment'] = met_dict['compartment']
        #args['charge'] = met_dict['charge']
        args['reactions'] = reactions


    def add_gene_data(self, args, gene_id):
        gene_dict = self.simulator.get_gene(gene_id)

        reactions = {rxn: self.mewpy_rxn_instance(rxn) for rxn in gene_dict['reactions']}

        args['types'] = set(['gene'])
        args['identifier'] = gene_id
        args['name'] = gene_dict['name']
        args['aliases'] = {gene_id, gene_dict['name']}
        args['model'] = self
        args['reactions'] = reactions


    def cobra_met_to_mewpy_met(self, met_id):
        args = {}
        self.add_metabolite_data(args, met_id)
        self.add_reg_data(args, met_id)

        return Variable.from_types(**args)


    def cobra_rxn_to_mewpy_rxn(self, rxn_id):
        args = {}
        self.add_reaction_data(args, rxn_id)
        self.add_reg_data(args, rxn_id)

        return Variable.from_types(**args)
    
    def cobra_gene_to_mewpy_gene(self, gene_id):
        args = {}
        self.add_gene_data(args, gene_id)
        self.add_reg_data(args, gene_id)

        return Variable.from_types(**args)
    

    def mewpy_rxn_instance(self, rxn_id):
        args = {'types': {'reaction'}, 'identifier': rxn_id, 'model': self}
        return Variable.from_types(**args)

    def cobra_objective_to_objective(self):
        res = {}
        for rxn_id, coeff in self.simulator.objective.items():
            rxn = self.cobra_rxn_to_mewpy_rxn(rxn_id)
            res[rxn] = float(coeff)
        return res
    


    def get(self, identifier: Any, default=None) -> Union['Gene', 'Metabolite', 'Reaction']:
        
        if not self.parsed:

            if identifier in self._metabolites:
                return self._metabolites[identifier]

            elif identifier in self._reactions:
                return self._reactions[identifier]

            elif identifier in self._genes:
                return self._genes[identifier]
            
            else:
                return super(MetabolicModelWrapper, self).get(identifier=identifier, default=default)
        
        else:
        
            if identifier in self.simulator.reactions:
                return self.cobra_rxn_to_mewpy_rxn(identifier)
                    
            elif identifier in self.simulator.genes:
                return self.cobra_gene_to_mewpy_gene(identifier)

            elif identifier in self.simulator.metabolites:
                return self.cobra_met_to_mewpy_met(identifier)

            else:
                return super(MetabolicModelWrapper, self).get(identifier=identifier, default=default) 


    def clean(self):
        self._metabolites = {}
        self._genes = {}
        self._reactions = {}
        self._parsed = 1


    def clean_history(self):
        super(MetabolicModelWrapper, self).clean_history()
        self.clean()



    def add(self,
            *variables: Union['Gene', 'Metabolite', 'Reaction'],
            comprehensive: bool = True,
            history: bool = True):
        
        if not self.parsed:
            for variable in variables:
                if 'gene' in variable.types:
                    self._genes[variable.id] = variable

                if 'metabolite' in variable.types:
                    self._metabolites[variable.id] = variable

                if 'reaction' in variable.types:
                    if comprehensive:

                        for metabolite in variable.yield_metabolites():
                            self._metabolites[metabolite.id] = metabolite

                        for gene in variable.yield_genes():
                            self._genes[gene.id] = gene

                    self._reactions[variable.id] = variable

        else:
            for variable in variables:
                if 'gene' in variable.types:
                    self.simulator.add_gene(variable.id, variable.name)

                if 'metabolite' in variable.types:
                    self.simulator.add_metabolite(variable.id,
                                                formula=variable.formula,
                                                name=variable.name,
                                                compartment=variable.compartment)

                if 'reaction' in variable.types:
                    self.simulator.add_reaction(id=variable.id,
                                                name=variable.name,
                                                stoichiometry=variable.stoichiometry,
                                                lb=variable.lower_bound,
                                                ub=variable.upper_bound,
                                                gpr=variable.gene_protein_reaction_rule)

        return super(MetabolicModelWrapper, self).add(*variables, comprehensive=comprehensive, history=history)
    

    def remove(self,
               *variables: Union['Gene', 'Metabolite', 'Reaction'],
               remove_orphans: bool = False,
               history: bool = True):
        for variable in variables:

            if 'gene' in variable.types:
                pass

            if 'metabolite' in variable.types:
                pass

            if 'reaction' in variable.types:
                self.simulator.remove_reaction(variable)
        
        return super(MetabolicModelWrapper, self).remove(*variables, remove_orphans=remove_orphans, history=history)


    def wrapper_simulation(self,  method='fba') -> Solution:
        if method == 'fba':
            res = self.simulator.simulate()

        elif method == 'pfba':
            res = self.simulator.simulate(method='pFBA')

        values = dict(res.fluxes)
        status = res.status
        fobj=res.objective_value
        message = status.value

        return Solution(fobj=fobj,values=values,status=status, message=message)


    def wrapper_fva(self,
                    fraction=0.9,
                    reactions: Sequence[str] = None,
                    objective=None,
                    constraints=None):
        if objective != None:
            raise ValueError("Cobrapy does not support changing the objective in FVA simulations")

        res = self.simulator.FVA(reactions=reactions, obj_frac=fraction, constraints=constraints, format='df')
        return res.rename(columns ={'Minimum': 'minimum', "Maximum":'maximum'})
        
