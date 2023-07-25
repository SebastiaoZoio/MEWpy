from typing import Union, TYPE_CHECKING

from .engine import Engine

from mewpy.io.dto import DataTransferObject
from mewpy.germ.variables.variable import Variable


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

            return CobraModelWrapper(identifier=identifier, model=self.io)

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


        for rxn in self.dto.cobra_model.reactions:
            self.variables[rxn.id].add('reaction')

        for met in self.dto.cobra_model.metabolites:
                self.variables[met.id].add('metabolite')

        for gene in self.dto.cobra_model.genes:
                self.variables[gene.id].add('gene')
    
    
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

        for var_id, types in variables.items():
            if len(types) > 1:
                if 'reaction' in types:
                    args = {}
                    model.add_reaction_data(args, var_id)
                    args['types'] = types
                    rxn = Variable.from_types(**args)

                    model.add_init_var(rxn)

                elif 'metabolite' in types:
                    args = {}
                    model.add_metabolite_data(args, var_id)
                    args['types'] = types
                    met = Variable.from_types(**args)

                    model.add_init_var(met)

                elif 'gene' in types:
                    args = {}
                    model.add_gene_data(args, var_id)
                    args['types'] = types
                    gene = Variable.from_types(**args)

                    model.add_init_var(gene)

        model.initializing = 1     

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
