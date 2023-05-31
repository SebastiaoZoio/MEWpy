from typing import Any, Sequence, TYPE_CHECKING, Generator, Union, List
from .metabolic_wrapper import MetabolicModelWrapper

from mewpy.germ.models.serialization import serialize
from mewpy.util.utilities import generator
from mewpy.util.history import HistoryManager

from ...solvers.solution import Solution

from cobra.core.solution import Solution as Cobra_Solution
from cobra.flux_analysis import pfba, flux_variability_analysis
from enum import Enum

from src.mewpy.germ.variables import Metabolite, Gene, Reaction
from cobra.core.reaction import Reaction as CobraRxn
from cobra.core.metabolite import Metabolite as CobraMet
from cobra.core.gene import Gene as CobraGene
from cobra.core import GPR as CobraGPR
from cobra.manipulation.delete import prune_unused_reactions, prune_unused_metabolites, remove_genes

from mewpy.germ.variables.variable import Variable
from mewpy.simulation.cobra import Simulation

from src.mewpy.io.engines.engines_utils import build_symbolic
from mewpy.germ.algebra import Expression
from mewpy.util.utilities import AttrDict


if TYPE_CHECKING:
    from mewpy.germ.algebra import Expression
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
            self._metabolites = {cobra_met.id: self.cobra_met_to_mewpy_met(cobra_met) for cobra_met in self.met_model.metabolites}
        return self._metabolites
    
    @property
    def genes(self):
        if not self._genes:
            self._genes = {cobra_gene.id: self.cobra_gene_to_mewpy_gene(cobra_gene) for cobra_gene in self.met_model.genes}
        return self._genes

    @property
    def reactions(self):
        if not self._reactions:
            self._reactions = {cobra_rxn.id: self.cobra_rxn_to_mewpy_rxn(cobra_rxn) for cobra_rxn in self.met_model.reactions}
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
                reg = self.regulators[id]
                add_regulator_data(args, reg)
    

            if id in self.targets:
                tar = self.targets[id]
                add_target_data(args, tar)



    def add_reaction_data(self, args:dict, cobra_rxn: CobraRxn):

        symbolic, warning = build_symbolic(expression=cobra_rxn.gene_reaction_rule)
        genes = {gene.id: self.cobra_gene_to_mewpy_gene(gene) for gene in cobra_rxn.genes}
        gpr = Expression(symbolic=symbolic, variables=genes)

        stoichiometry = {self.cobra_met_to_mewpy_met(met): cobra_rxn.metabolites[met] for met in cobra_rxn.metabolites}
        
        args['types']=set(['reaction'])
        args['identifier']=cobra_rxn.id
        args['name']=cobra_rxn.name
        args['aliases']={cobra_rxn.id, cobra_rxn.name}
        args['model']=self
        args['bounds']=cobra_rxn.bounds
        args['stoichiometry']=stoichiometry
        args['gpr']=gpr


    
    def add_metabolite_data(self, args:dict, cobra_met: CobraMet):

        reactions = {rxn.id: mewpy_rxn_instance(rxn, self) for rxn in cobra_met.reactions}

        args['types'] = set(['metabolite'])
        args['identifier'] = cobra_met.id
        args['name'] = cobra_met.name
        args['aliases'] = {cobra_met.id, cobra_met.name}
        args['model'] = self
        args['formula'] = cobra_met.formula
        args['charge'] = cobra_met.charge
        args['compartment'] = cobra_met.compartment
        args['reactions'] = reactions
        
        

    def add_gene_data(self, args, cobra_gene: CobraGene):

        reactions = {rxn.id: mewpy_rxn_instance(rxn, self) for rxn in cobra_gene.reactions}

        args['types'] = set(['gene'])
        args['identifier'] = cobra_gene.id
        args['name'] = cobra_gene.name
        args['aliases'] = {cobra_gene.id, cobra_gene.name}
        args['model'] = self
        args['reactions'] = reactions


    def cobra_met_to_mewpy_met(self, cobra_met: CobraMet):
        args = {}
        self.add_metabolite_data(args, cobra_met)
        self.add_reg_data(args, cobra_met.id)

        return Variable.from_types(**args)

    def cobra_rxn_to_mewpy_rxn(self, cobra_rxn: AttrDict):
        args = {}
        self.add_reaction_data(args, cobra_rxn)
        self.add_reg_data(args, cobra_rxn.id)

        return Variable.from_types(**args)
    
    def cobra_gene_to_mewpy_gene(self, cobra_gene: CobraGene):
        args = {}
        self.add_gene_data(args, cobra_gene)
        self.add_reg_data(args, cobra_gene.id)

        return Variable.from_types(**args)
    
    def cobra_objective_to_objective(self):
        res = {}
        for arg in self.met_model.objective.expression.args:
            coef, reaction_id = str(arg).split('*')
            try:
                cobra_rxn = self.met_model.reactions.get_by_id(reaction_id)
                reaction = self.cobra_rxn_to_mewpy_rxn(cobra_rxn)
                res[reaction] = float(coef)
            except:
                continue

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
        
            if identifier in self.met_model.reactions:
                return self.cobra_rxn_to_mewpy_rxn(self.met_model.reactions.get_by_id(identifier))
                    
            elif identifier in self.met_model.genes:
                return self.cobra_gene_to_mewpy_gene(self.met_model.genes.get_by_id(identifier))

            elif identifier in self.met_model.metabolites:
                return self.cobra_met_to_mewpy_met(self.met_model.metabolites.get_by_id(identifier))

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
                remove_genes(self.met_model, [variable])

            if 'metabolite' in variable.types:
                self.met_model.remove_metabolites([variable])

            if 'reaction' in variable.types:
                self.met_model.remove_reactions([variable])
            
            if remove_orphans:
                prune_unused_metabolites(self.met_model)
                prune_unused_reactions(self.met_model)

        
        return super(MetabolicModelWrapper, self).remove(*variables, remove_orphans=remove_orphans, history=history)



    def wrapper_simulation(self,  method='fba') -> Solution:
        if method == 'fba':
            cobra_solution = self.met_model.optimize()
            solution = cobra_solution_to_solution(cobra_solution)

            return solution

        elif method == 'pfba':
            cobra_solution = pfba(self.met_model)
            solution = cobra_solution_to_solution(cobra_solution)

            return solution


    def wrapper_fva(self,
                    fraction=1.0,
                    reactions: Sequence[str] = None,
                    objective=None,
                    constraints=None):
        if objective != None:
            raise ValueError("Cobrapy does not support changing the objective in FVA simulations")
        if constraints != None:
            raise ValueError("Cobrapy does not support adding constraints in FVA simulations")
        
        return flux_variability_analysis(model=self.met_model, fraction_of_optimum=fraction, reaction_list=reactions)









## aux functions


class Status(Enum):
    """ Enumeration of possible solution status. """
    OPTIMAL = 'Optimal'
    UNKNOWN = 'Unknown'
    SUBOPTIMAL = 'Suboptimal'
    UNBOUNDED = 'Unbounded'
    INFEASIBLE = 'Infeasible'
    INF_OR_UNB = 'Infeasible or Unbounded'

def get_status(status: str) -> Status:
    if status == 'optimal':
        return Status.OPTIMAL
    elif status == 'suboptimal':
        return Status.SUBOPTIMAL
    elif status == 'unbounded':
        return Status.UNBOUNDED
    elif status == 'infeasible':
        return Status.INFEASIBLE
    elif status == 'infeasible_or_unbounded':
        return Status.INF_OR_UNB
    else:
        return Status.UNKNOWN
    

def cobra_solution_to_solution(wrapper_solution: Cobra_Solution):
    fobj = wrapper_solution.objective_value
    message = wrapper_solution.status
    values = wrapper_solution.fluxes.to_dict()
    status = get_status(wrapper_solution.status)

    return Solution(status=status, message=message, fobj=fobj, values=values)


def mewpy_rxn_instance(cobra_rxn:CobraRxn, model):
    rxn = Reaction(identifier= cobra_rxn.id)
    rxn.update(name= cobra_rxn.name, aliases={cobra_rxn.id, cobra_rxn.name}, model=model)

    return rxn
  

def add_regulator_data(args:dict, regulator):
    args['types'].add('regulator')
    args['interactions'] = regulator.interactions


def add_target_data(args:dict, target):
    args['types'].add('target')
    args['interaction'] = target.interaction



def mewpy_gene_to_cobra_gene(gene : Gene) -> CobraGene:
    return CobraGene(gene.id, name=gene.name, functional=gene.is_active)


def mewpy_met_to_cobra_met(met : Metabolite) -> CobraMet:    
    return CobraMet(met.id,formula=met.formula,name=met.name,compartment=met.compartment, charge=met.charge)


def mewpy_rxn_to_cobra_rxn(rxn : Reaction) -> CobraRxn:
    cobra_rxn = CobraRxn(rxn.id, name=rxn.name, lower_bound=rxn.lower_bound, upper_bound=rxn.upper_bound)

    cobra_rxn.gpr = CobraGPR.from_string(rxn.gpr.to_string())    

    mets = {}
    for met in rxn.metabolites.values():
        #get coefficient
        mets[mewpy_met_to_cobra_met(met)] = rxn.stoichiometry[met]
    cobra_rxn.add_metabolites(mets)
        
    return cobra_rxn