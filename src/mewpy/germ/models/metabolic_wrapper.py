from typing import TYPE_CHECKING, Any, Union, Generator, List
from .model import Model
from mewpy.germ.models.serialization import serialize
from mewpy.util.utilities import generator
from mewpy.util.history import HistoryManager


if TYPE_CHECKING:
    from mewpy.germ.algebra import Expression
    from mewpy.germ.variables import Gene, Metabolite, Reaction

class MetabolicModelWrapper(Model, model_type='metabolic_wrapper', register=True, constructor=True, checker=True):
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
    def metabolites(self):
        pass

    @property
    def genes(self):
        pass
    

    
    @property
    def reactions(self):
        pass
    

    @property
    def objective(self):
        pass

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
    

    def wrapper_simulation(self, method='fba'):
        pass


    def get(self, identifier: Any, default=None) -> Union['Gene', 'Metabolite', 'Reaction']:
        return super(MetabolicModelWrapper, self).get(identifier=identifier, default=default)

    def wrapper_fva(self,
                    fraction,
                    reactions,
                    objective,
                    constraints):
        pass
