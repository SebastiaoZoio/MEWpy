
from .metabolic import MetabolicModel

from mewpy.simulation.reframed import Simulation
from mewpy.germ.variables.variable import Variable
from mewpy.io.engines.engines_utils import build_symbolic

from reframed import CBModel as ReframedModel

REFRAMED_METHODS = ('FBA', 'pFBA', 'FVA')

class ReframedModelWrapper(MetabolicModel, model_type='reframed_wrapper', register=True, constructor=True, checker=True):

    def __init__(self,
                 identifier,
                 **kwargs):

        super().__init__(identifier,
                         **kwargs)
        
    
    
    @property
    def simulator(self):
        if not self._simulator:
            reframed_model = ReframedModel("reframed model")
            self._simulator = Simulation(reframed_model)

        return self._simulator
    


    @simulator.setter
    def simulator(self, simulator):
        self._simulator = simulator

    def set_simulator(self, met_model):
        self.simulator = Simulation(met_model)
        
    
    def clean(self):
        self._metabolites = {}
        self._genes = {}
        self._reactions = {}

    def clean_history(self):
        super(MetabolicModel, self).clean_history()
        self.destroy_init_vars()

    def has_external_method(self, method:str):
        return method in REFRAMED_METHODS
    

    def build_variable(self, args):
        return Variable.from_types(**args)
    
    def build_symbolic(self, expression):
        return build_symbolic(expression)
