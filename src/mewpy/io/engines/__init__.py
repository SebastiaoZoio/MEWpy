from enum import Enum

from .boolean_csv import BooleanRegulatoryCSV
from . co_expression_csv import CoExpressionRegulatoryCSV
from .target_regulator_csv import TargetRegulatorRegulatoryCSV
from .json import JSON
from .metabolic_sbml import MetabolicSBML
from .regulatory_sbml import RegulatorySBML
from .cobra_model_engine import CobraModelEngine
from .reframed_model_engine import ReframedModelEngine


class Engines(Enum):
    """
    List of all engines to read/write files/models

    """

    BooleanRegulatoryCSV = BooleanRegulatoryCSV
    CoExpressionRegulatoryCSV = CoExpressionRegulatoryCSV
    TargetRegulatorRegulatoryCSV = TargetRegulatorRegulatoryCSV
    RegulatorySBML = RegulatorySBML
    MetabolicSBML = MetabolicSBML
    CobraModel = CobraModelEngine
    ReframedModel = ReframedModelEngine
    JSON = JSON
    CobraModelEngine = CobraModelEngine
    ReframedModelEngine = ReframedModelEngine

    @classmethod
    def has_engine(cls, engine):
        try:
            return cls[engine.__name__]

        except KeyError:

            return False

    @classmethod
    def get(cls, engine, default=None):
        try:
            return cls[engine.__name__]

        except KeyError:

            return default
