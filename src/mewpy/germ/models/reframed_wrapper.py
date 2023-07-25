
from .metabolic import MetabolicModel

from mewpy.germ.variables.variable import Variable
from mewpy.io.engines.engines_utils import build_symbolic


class ReframedModelWrapper(MetabolicModel, model_type='reframed_wrapper', register=True, constructor=True, checker=True):

    def __init__(self,
                 identifier,
                 **kwargs):

        super().__init__(identifier,
                         **kwargs)
        
        self.external_methods = ('FBA', 'pFBA', 'FVA')
  

    def build_variable(self, args):
        return Variable.from_types(**args)
    
    def build_symbolic(self, expression):
        return build_symbolic(expression)
