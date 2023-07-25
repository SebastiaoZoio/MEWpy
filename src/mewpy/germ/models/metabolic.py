from typing import TYPE_CHECKING, Any, Sequence, Union, Generator, List, Dict, Tuple, Set
from .model import Model
from mewpy.germ.models.serialization import serialize
from mewpy.util.utilities import generator
from mewpy.util.history import HistoryManager, recorder
from mewpy.solvers.solution import Solution, Status, SStatus
from mewpy.germ.algebra import Expression
from mewpy.simulation.simulator import get_simulator

if TYPE_CHECKING:
    
    from mewpy.germ.variables import Gene, Metabolite, Reaction

class MetabolicModel(Model, model_type='metabolic', register=True, constructor=True, checker=True):
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
                 model=None,
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
        # compartments attribute can be shared across the children, thus name mangling
        self.__compartments = {}
        self._genes = {}
        self._metabolites = {}
        self._objective = {}
        self._reactions = {}
        self.simulator = None

        # auxiliar variables to build integrated models
        self.initializing = 0
        self.init_rxns = {}
        self.init_mets = {}
        self.init_genes = {}

        super().__init__(identifier,
                         **kwargs)
        
        self.set_simulator(model)
        self.external_methods = ()


    @classmethod
    def create(cls, tool, id, **kwargs):
        if tool == "mewpy":
            from .mewpy_metabolic import MewpyMetabolicModel
            return MewpyMetabolicModel(id, **kwargs)
        elif tool == "cobra":
            from .cobra_wrapper import CobraModelWrapper
            return CobraModelWrapper(id, **kwargs)
        else:
            raise ValueError("Invalid tool specified")


    @classmethod
    def with_tool(cls, tool, id, **kwargs):
        return cls.create(tool, id, **kwargs)

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
        _types = {MetabolicModel.model_type}

        _types.update(super(MetabolicModel, self).types)

        return _types

    # -----------------------------------------------------------------------------
    # Static attributes
    # -----------------------------------------------------------------------------

    @property
    def metabolites(self) -> List[str]:
        if not self.simulator:
            return {}
        if not self._metabolites:
            self._metabolites = dict((m, self.external_met_to_mewpy_met(m)) for m in self.simulator.metabolites)
        return self._metabolites
        
    @property
    def genes(self) -> List[str]:
        if not self.simulator:
            return {}

        if not self._genes:
            self._genes = dict((g, self.external_gene_to_mewpy_gene(g)) for g in self.simulator.genes)
        return self._genes

    @property
    def reactions(self) -> List[str]:
        if not self.simulator:
            return {}

        if not self._reactions:
            self._reactions = dict((r, self.external_rxn_to_mewpy_rxn(r)) for r in self.simulator.reactions)
        return self._reactions
    

    @property
    def objective(self):
        return self.external_objective_to_objective()
    

    @property
    def exchanges(self):
        return dict((ext_rxn, self.external_rxn_to_mewpy_rxn(ext_rxn)) for ext_rxn in self.simulator.get_exchange_reactions())
    
    @property
    def demands(self):
        return dict((ext_rxn, self.external_rxn_to_mewpy_rxn(ext_rxn)) for ext_rxn in self.simulator.get_demand_reactions())

    @property
    def sinks(self):
        return dict((ext_rxn, self.external_rxn_to_mewpy_rxn(ext_rxn)) for ext_rxn in self.simulator.get_sink_reactions())
    

    @property
    def contexts(self) -> List[HistoryManager]:
        return super(MetabolicModel, self).contexts
    

    @property
    def compartments(self) -> Dict[str, str]:
        """
        It returns a dictionary with the compartments of the model. The key is the compartment identifier a
        nd the value is the compartment name.
        To retrieve an iterator with the compartments use `yield_compartments` method.
        Note that the compartments attribute retrieves a copy of the compartments' container.
        To update the compartments container set new `compartments`.
        :return:
        """

        compartments = {met.compartment: self.__compartments.get(met.compartment, '')
                        for met in self.yield_metabolites()
                        if met.compartment is not None}

        compartments.update(self.__compartments)

        compartments.update(super(MetabolicModel, self).compartments)

        return compartments

    # -----------------------------------------------------------------------------
    # Static attributes setters
    # -----------------------------------------------------------------------------
    def set_simulator(self, model):
        if model:
            self.simulator = get_simulator(model)


    @compartments.setter
    @recorder
    def compartments(self, value: Dict[str, str]):
        """
        It sets the compartments of the model. The key is the compartment identifier a
        nd the value is the compartment name.
        :param value: a dictionary with the compartments of the model
        :return:
        """
        if not value:
            value = {}

        self.__compartments.update(value)

    
    @genes.setter
    @recorder
    def genes(self, value: Dict[str, 'Gene']):
        """
        It sets the genes of the model. The key is the gene identifier and the value is the `Gene` object.
        :param value: a dictionary with the genes of the model
        :return:
        """
        if not value:
            value = {}

        self.add(*value.values(), history=False)

    @metabolites.setter
    @recorder
    def metabolites(self, value: Dict[str, 'Metabolite']):
        """
        It sets the metabolites of the model. The key is the metabolite identifier and the value is the `Metabolite`
        object.
        :param value: a dictionary with the metabolites of the model
        :return:
        """
        if not value:
            value = {}

        self.add(*value.values(), history=False)

    
    @reactions.setter
    @recorder
    def reactions(self, value: Dict[str, 'Reaction']):
        """
        It sets the reactions of the model. The key is the reaction identifier and the value is the `Reaction` object.
        :param value: a dictionary with the reactions of the model
        :return:
        """
        if not value:
            value = {}

        self.add(*value.values(), history=False)


    @objective.setter
    @recorder
    def objective(self, objective):
        self.simulator.objective = objective


    # -----------------------------------------------------------------------------
    # Generators
    # -----------------------------------------------------------------------------

    def yield_genes(self):
        """
        It yields the genes of the model.
        :return: a generator with the genes of the model
        """
        #for g in self.genes:
        #    yield self.external_gene_to_mewpy_gene(g)

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
        #for m in self.metabolites:
        #    yield self.external_met_to_mewpy_met(m)

        return generator(self.metabolites)

    def yield_reactions(self) -> Generator['Reaction', None, None]:
        """
        It yields the reactions of the model.
        :return: a generator with the reactions of the model
        """
        #for r in self.reactions:
        #    yield self.external_rxn_to_mewpy_rxn(r)

        return generator(self.reactions)
    


    # -----------------------------------------------------------------------------
    # Operations/Manipulations
    # -----------------------------------------------------------------------------

    def get(self, identifier: Any, default=None) -> Union['Gene', 'Metabolite', 'Reaction']:

        if self.initializing:
            return self.get_init_var(identifier)

        if identifier in self.simulator.reactions:
            return self.external_rxn_to_mewpy_rxn(identifier)
                    
        elif identifier in self.simulator.genes:
            return self.external_gene_to_mewpy_gene(identifier)

        elif identifier in self.simulator.metabolites:
            return self.external_met_to_mewpy_met(identifier)

        else:
            return super(MetabolicModel, self).get(identifier=identifier, default=default) 
        


    def add(self,
            *variables: Union['Gene', 'Metabolite', 'Reaction'],
            comprehensive: bool = True,
            history: bool = True):
        
        self.clean()

        for variable in variables:
            if 'gene' in variable.types:
                if variable.id not in self.simulator.genes:
                    self.simulator.add_gene(variable.id, variable.name)
            if 'metabolite' in variable.types:
                if variable.id not in self.simulator.metabolites:
                    self.simulator.add_metabolite(variable.id,
                                            formula=variable.formula,
                                            name=variable.name,
                                            compartment=variable.compartment)
            if 'reaction' in variable.types:
                if variable.id not in self.simulator.reactions:
                    self.simulator.add_reaction(rxn_id=variable.id,
                                            name=variable.name,
                                            stoichiometry={r.id: v for r , v in variable.stoichiometry.items()},
                                            lb=variable.lower_bound,
                                            ub=variable.upper_bound,
                                            gpr=variable.gene_protein_reaction_rule)


        return super(MetabolicModel, self).add(*variables, comprehensive=comprehensive, history=history)
    

    def remove(self,
               *variables: Union['Gene', 'Metabolite', 'Reaction'],
               remove_orphans: bool = False,
               history: bool = True):
        
        self.clean()

        for variable in variables:
            if 'reaction' in variable.types:
                self.simulator.remove_reaction(variable.id)

        
        return super(MetabolicModel, self).remove(*variables, remove_orphans=remove_orphans, history=history)


    def update(self,
               compartments: Dict[str, str] = None,
               objective: Dict['Reaction', Union[float, int]] = None,
               variables: Union[List[Union['Gene', 'Metabolite', 'Reaction']],
                                Tuple[Union['Gene', 'Metabolite', 'Reaction']],
                                Set[Union['Gene', 'Metabolite', 'Reaction']]] = None,
               **kwargs):
        """
        It updates the model with relevant information, namely the compartments, objective and variables.

        :param compartments: the compartments to be updated
        :param objective: the objective to be updated
        :param variables: the variables to be updated
        :param kwargs: additional arguments
        :return:
        """

        self.clean()

        if compartments is not None:
            self.compartments = compartments

        if variables is not None:
            self.add(*variables)

        if objective is not None:
            self.objective = objective

        super(MetabolicModel, self).update(**kwargs)
    
    def clean(self):
        self._metabolites = {}
        self._genes = {}
        self._reactions = {}


    def clean_history(self):
        super(MetabolicModel, self).clean_history()
        self.destroy_init_vars()



    def get_init_var(self, identifier: Any) -> Union['Gene', 'Metabolite', 'Reaction']:
        if identifier in self.init_mets:
            return self.init_mets[identifier]

        elif identifier in self.init_rxns:
            return self.init_rxns[identifier]

        elif identifier in self.init_genes:
            return self.init_genes[identifier]


    def add_init_var(self, variable: Union['Gene', 'Metabolite', 'Reaction']):
        if 'gene' in variable.types:
            self.init_genes[variable.id] = variable

        if 'metabolite' in variable.types:
            self.init_mets[variable.id] = variable

        if 'reaction' in variable.types:
            self.init_rxns[variable.id] = variable


    def destroy_init_vars(self):
        self.init_genes = {}
        self.init_mets = {}
        self.init_rxns = {}
        self.initializing = 0

    # -----------------------------------------------------------------------------
    # Variable Converters
    # -----------------------------------------------------------------------------
    
    def external_met_to_mewpy_met(self, met_id):
        args = {}
        self.add_metabolite_data(args, met_id)
        self.add_reg_data(args, met_id)

        return self.build_variable(args)


    def external_rxn_to_mewpy_rxn(self, rxn_id):
        args = {}
        self.add_reaction_data(args, rxn_id)
        self.add_reg_data(args, rxn_id)

        return self.build_variable(args)
    
    def external_gene_to_mewpy_gene(self, gene_id):
        args = {}
        self.add_gene_data(args, gene_id)
        self.add_reg_data(args, gene_id)

        return self.build_variable(args)
    
    def external_objective_to_objective(self):
        res = {}
        for rxn_id, coeff in self.simulator.objective.items():
            rxn = self.external_rxn_to_mewpy_rxn(rxn_id)
            res[rxn] = float(coeff)
        return res
    

    def mewpy_rxn_instance(self, rxn_id):
        args = {'types': {'reaction'}, 'identifier': rxn_id, 'model': self}
        return self.build_variable(args)


    def mewpy_met_instance(self, met_id):
        args = {'types': {'metabolite'}, 'identifier': met_id, 'model': self}
        met_data = self.simulator.get_metabolite(met_id)
        args['compartment'] = met_data['compartment']
        args['formula'] = met_data['formula']
        return self.build_variable(args)
    

    def mewpy_gene_instance(self, gene_id):
        args = {'types': {'gene'}, 'identifier': gene_id, 'model': self}
        return self.build_variable(args)



    def add_reaction_data(self, args:dict, rxn_id):
        rxn_dict = self.simulator.get_reaction(rxn_id)

        symbolic, warning = self.build_symbolic(expression=rxn_dict['gpr'])
        genes = {gene.name: self.mewpy_gene_instance(gene.name) for gene in symbolic.atoms(symbols_only=True)}
        gpr = Expression(symbolic=symbolic, variables=genes)

        stoichiometry = {self.mewpy_met_instance(met_id): coeff for met_id, coeff in rxn_dict['stoichiometry'].items()}

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

        reactions = {rxn: self.external_rxn_to_mewpy_rxn(rxn) for rxn in self.simulator.get_metabolite_reactions(met_id)}

        args['types'] = set(['metabolite'])
        args['identifier'] = met_id
        args['name'] = met_dict['name']
        args['aliases'] = {met_id, met_dict['name']}
        args['model'] = self
        args['formula'] = met_dict['formula']
        args['compartment'] = met_dict['compartment']
        args['charge'] = met_dict['charge']
        args['reactions'] = reactions


    def add_gene_data(self, args, gene_id):
        gene_dict = self.simulator.get_gene(gene_id)

        reactions = {rxn: self.external_rxn_to_mewpy_rxn(rxn) for rxn in gene_dict['reactions']}

        args['types'] = set(['gene'])
        args['identifier'] = gene_id
        args['name'] = gene_dict['name']
        args['aliases'] = {gene_id, gene_dict['name']}
        args['model'] = self
        args['reactions'] = reactions


    def add_reg_data(self, args, id):
        if self.is_regulatory():
            if id in self.regulators:
                args['types'].add('regulator')
                args['interactions'] = self.regulators[id].interactions
        
            if id in self.targets:
                args['types'].add('target')
                args['interaction'] = self.targets[id].interaction
    


    # -----------------------------------------------------------------------------
    # Simulation Converters
    # -----------------------------------------------------------------------------    

    def wrapper_simulation(self,  method='fba', constraints=None, slim=False) -> Solution:
        if method == 'fba':
            res = self.simulator.simulate(constraints=constraints)

        elif method == 'pfba':
            res = self.simulator.simulate(method='pFBA', constraints=constraints)
        
        if slim:
            return res.objective_value

        values=None
        if res.fluxes:
            values = dict(res.fluxes)
        status = status_mapping[res.status]
        fobj=res.objective_value
        message = status.value

        return Solution(fobj=fobj,values=values,status=status, message=message)


    def wrapper_fva(self,
                    fraction=0.9,
                    reactions: Sequence[str] = None,
                    constraints=None):
    
        res = self.simulator.FVA(reactions=reactions, obj_frac=fraction, constraints=constraints, format='df')
        return res.rename(columns ={'Minimum': 'minimum', "Maximum":'maximum'})

    def single_gene_deletion(self, genes: Union[List[str], None] = None):
        return self.simulator.gene_deletion(genes=genes)
    
    def single_reaction_deletion(self, reactions: Union[List[str], None] = None):
        return self.simulator.reaction_deletion(reactions=reactions)

    def has_external_method(self, method:str):
        return method in self.external_methods
    


status_mapping = {
    SStatus.OPTIMAL: Status.OPTIMAL,
    SStatus.UNKNOWN: Status.UNKNOWN,
    SStatus.SUBOPTIMAL: Status.SUBOPTIMAL,
    SStatus.UNBOUNDED: Status.UNBOUNDED,
    SStatus.INFEASIBLE: Status.INFEASIBLE,
    SStatus.INF_OR_UNB: Status.INF_OR_UNB
}