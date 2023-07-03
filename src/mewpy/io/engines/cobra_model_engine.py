from functools import partial
from typing import Union, TYPE_CHECKING

from .engine import Engine

from mewpy.io.dto import VariableRecord, DataTransferObject, CompartmentRecord, FunctionTerm
from .engines_utils import build_symbolic, expression_warning, cobra_warning


from mewpy.germ.algebra import Expression
from ...germ.models.cobra_wrapper import CobraModelWrapper

if TYPE_CHECKING:
    from ...germ.models import RegulatoryModel, Model, MetabolicModel


class CobraModelEngine(Engine):
    def __init__(self, io, config, model=None):
        """
        Engine for COBRApy constraint-based metabolic models
        """
        super().__init__(io, config, model)

    @property
    def model_type(self):
        return 'cobra_wrapper'


    @property
    def model(self):

        if self._model is None:
            identifier = self.get_identifier()

            return CobraModelWrapper(identifier=identifier)

        return self._model

    def open(self, mode='r'):

        self._dto = DataTransferObject()

        if not hasattr(self.io, 'reactions'):
            raise OSError(f'{self.io} is not a valid input. Provide a cobrapy model')

        self.dto.cobra_model = self.io

        self.dto.id = self.get_identifier()

        self.dto.name = self.dto.cobra_model.name

    def parse(self):

        if self.dto is None:
            raise OSError('Model is not open')

        if self.dto.id is None:
            raise OSError('Model is not open')

        if self.dto.cobra_model is None:
            raise OSError('Model is not open')

        # -----------------------------------------------------------------------------
        # Compartments
        # -----------------------------------------------------------------------------

        self.dto.compartments = {c_id: CompartmentRecord(id=c_id, name=c_name)
                                 for c_id, c_name in self.dto.cobra_model.compartments.items()}

        # -----------------------------------------------------------------------------
        # Reactions
        # -----------------------------------------------------------------------------
        processed_metabolites = set()
        processed_genes = set()

        for rxn in self.dto.cobra_model.reactions:

            # -----------------------------------------------------------------------------
            # Genes
            # -----------------------------------------------------------------------------

            symbolic, warning = build_symbolic(expression=rxn.gene_reaction_rule)

            if warning:
                self.warnings.append(partial(expression_warning, warning))

            genes = {}

            for symbol in symbolic.atoms(symbols_only=True):
                self.variables[symbol.name].add('gene')

                gene_record = VariableRecord(id=symbol.name,
                                             name=symbol.name,
                                             aliases={symbol.name, symbol.value})

                genes[symbol.name] = gene_record

                processed_genes.add(symbol.name)

                self.dto.genes[symbol.name] = gene_record

            # -----------------------------------------------------------------------------
            # Metabolites
            # -----------------------------------------------------------------------------

            stoichiometry = {}
            metabolites = {}

            for met in rxn.metabolites:
                met_record = VariableRecord(id=met.id,
                                            name=met.name,
                                            aliases={met.id, met.name},
                                            compartment=met.compartment,
                                            charge=met.charge,
                                            formula=met.formula)

                metabolites[met.id] = met_record

                coef = rxn.metabolites[met]

                stoichiometry[met.id] = coef

                processed_metabolites.add(met.id)

                self.variables[met.id].add('metabolite')

                self.dto.metabolites[met.id] = met_record

            # -----------------------------------------------------------------------------
            # GPR Function term
            # -----------------------------------------------------------------------------
            function_term = FunctionTerm(id='gpr_term', symbolic=symbolic, coefficient=1)

            # -----------------------------------------------------------------------------
            # Reaction
            # -----------------------------------------------------------------------------
            reaction_record = VariableRecord(id=rxn.id,
                                             name=rxn.name,
                                             aliases={rxn.id, rxn.name},
                                             bounds=rxn.bounds,
                                             genes=genes,
                                             gpr=function_term,
                                             stoichiometry=stoichiometry,
                                             metabolites=metabolites)

            self.variables[rxn.id].add('reaction')

            self.dto.reactions[rxn.id] = reaction_record

        for met in self.dto.cobra_model.metabolites:

            if met.id not in processed_metabolites:
                met_record = VariableRecord(id=met.id,
                                            name=met.name,
                                            aliases={met.id, met.name},
                                            compartment=met.compartment,
                                            charge=met.charge,
                                            formula=met.formula)

                self.variables[met.id].add('metabolite')

                self.dto.metabolites[met.id] = met_record

        for gene in self.dto.cobra_model.genes:

            if gene.id not in processed_genes:
                gene_record = VariableRecord(id=gene.id,
                                             name=gene.id,
                                             aliases={gene.id})

                self.variables[gene.id].add('gene')

                self.dto.genes[gene.id] = gene_record
    
    
    def read(self,
            model: Union['Model', 'MetabolicModel', 'RegulatoryModel'] = None,
            variables = None):
        
        if not model:
            model: Union['Model', 'MetabolicModel', 'RegulatoryModel'] = self.model

        if not variables:
            variables = self.variables

        if self.dto.id:
            model._id = self.dto.id

        if self.dto.name:
            model.name = self.dto.name

        model.set_simulator(self.dto.cobra_model)

        processed_metabolites = set()
        processed_genes = set()

        germ_reactions, germ_metabolites, germ_genes = get_germ_vars(variables, model, self.dto.cobra_model)


        for rxn_id in germ_reactions:
            rxn_record = self.dto.reactions[rxn_id]

            genes = {}

            for gene_id, gene_record in rxn_record.genes.items():

                gene, warning = gene_record.to_variable(model=model,
                                                        types=variables.get(gene_id, {'gene'}),
                                                        name=gene_record.name,
                                                        aliases=gene_record.aliases)

                if warning:
                    self.warnings.append(partial(cobra_warning, warning))

                genes[gene_id] = gene

                processed_genes.add(gene_id)

            stoichiometry = {}

            for met_id, met_record in rxn_record.metabolites.items():

                met, warning = met_record.to_variable(model=model,
                                                      types=variables.get(met_id, {'metabolite'}),
                                                      name=met_record.name,
                                                      aliases=met_record.aliases,
                                                      compartment=met_record.compartment,
                                                      charge=met_record.charge,
                                                      formula=met_record.formula)

                if warning:
                    self.warnings.append(partial(cobra_warning, warning))

                coef = rxn_record.stoichiometry[met_id]

                stoichiometry[met] = coef

                processed_metabolites.add(met_id)

            gpr = Expression(symbolic=rxn_record.gpr.symbolic, variables=genes)

            rxn, warning = rxn_record.to_variable(model=model,
                                                  types=variables.get(rxn_id, {'reaction'}),
                                                  bounds=rxn_record.bounds,
                                                  gpr=gpr,
                                                  stoichiometry=stoichiometry)

            if warning:
                self.warnings.append(partial(cobra_warning, warning))

        for met_id in germ_metabolites:
            met_record = self.dto.metabolites[met_id]

            if met_id not in processed_metabolites:
                met, warning = met_record.to_variable(model=model,
                                                      types=variables.get(met_id, {'metabolite'}),
                                                      name=met_record.name,
                                                      aliases=met_record.aliases,
                                                      compartment=met_record.compartment,
                                                      charge=met_record.charge,
                                                      formula=met_record.formula)

                if warning:
                    self.warnings.append(partial(cobra_warning, warning))


        for gene_id in germ_genes:
            gene_record = self.dto.genes[gene_id]

            if gene_id not in processed_genes:
                gene, warning = gene_record.to_variable(model=model,
                                                        types=variables.get(gene_id, {'gene'}),
                                                        name=gene_record.name,
                                                        aliases=gene_record.aliases)

                if warning:
                    self.warnings.append(partial(cobra_warning, warning))

        return model


    def write(self):
        pass

    def close(self):
        pass

    def clean(self):
        self._dto = None

    def get_identifier(self):

        if self.dto.cobra_model:
            return self.dto.cobra_model.id

        return 'model'



def get_germ_vars(variables, model, cobra_model):
    germ_sets = {key: set() for key in ['reaction', 'metabolite', 'gene']}

    for var_id, types in variables.items():
        if any(t in germ_sets for t in types) and len(types)>1:
            for t in types:
                if t in germ_sets:
                    germ_sets[t].add(var_id)
            model.add_reg_var(var_id, types)
    
    add_germ_mets_rxns(germ_sets['metabolite'], germ_sets['reaction'], cobra_model)
    add_germ_genes_rxns(germ_sets['gene'], germ_sets['reaction'], cobra_model)

    return germ_sets['reaction'], germ_sets['metabolite'], germ_sets['gene']


def add_germ_mets_rxns(germ_metabolites, germ_reactions, cobra_model):
    """Add reactions in which the regulatory metabolite is envolved"""
    for met in germ_metabolites:
        cobra_met = cobra_model.metabolites.get_by_id(met)
        for rxn in cobra_met.reactions:
            germ_reactions.add(rxn.id)


def add_germ_genes_rxns(germ_genes, germ_reactions, cobra_model):
    """Add reactions in which the regulatory gene is envolved"""
    for gene in germ_genes:
        cobra_gene = cobra_model.genes.get_by_id(gene)
        for rxn in cobra_gene.reactions:
            germ_reactions.add(rxn.id)