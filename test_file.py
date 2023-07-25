import warnings
warnings.filterwarnings("ignore", message="scipy._lib.messagestream.MessageStream size changed")
from termcolor import colored

from src.mewpy.io import Reader, Engines, read_model
from mewpy.germ.analysis import FBA, fva, SRFBA, RFBA
from cobra.io import read_sbml_model

FAILED = 0
SUCCESS = 1
DECIMAL_DIGITS = 6

e_coli = 1
yeast = 0

output_value = SUCCESS

def set_output_value(value):
    global output_value
    output_value = value

def round_(number):
    return round(number, DECIMAL_DIGITS)


def compare_dicts(dict1, dict2):
    """
    Compare the values of two dictionaries with the same keys.
    Prints a message for each key where the values are different.
    """
    for key in dict1.keys():
        if dict1[key] != dict2[key]:
            print(f"Values for key '{key}' are different:")
            print(f"  fba: {dict1[key]}")
            print(f"  new_fba: {dict2[key]}")
            

def compare_fva_dataframes(new_fva, fva):
    equal = True

    for (cobra_index, new_row), (mewpy_index, original_row) in zip(new_fva.iterrows(), fva.iterrows()):
        if (round_(new_row["minimum"]) != round_(original_row["minimum"]) or
            round_(new_row["maximum"]) != round_(original_row["maximum"])):
            print(colored(f"FVA Values for '{cobra_index}' differ", "red"))
            print(colored(f"  new fva: minimum->{new_row['minimum']} | maximum->{new_row['maximum']}", "blue"))
            print(colored(f"  orig fva: minimum->{original_row['minimum']} | maximum->{original_row['maximum']}", "blue"))
            equal = False
    return equal



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
    cobra_model = read_sbml_model(met_file)
    cobra_reader = Reader(Engines.CobraModelEngine, cobra_model)

    if reg_file:
        trn_reader = get_regulatory_reader(reg_file)
        return read_model(cobra_reader, trn_reader)

    return read_model(cobra_reader)



###########################################################################


def do_fba_asserts(original_model, new_model):

    fba = FBA(original_model).build().optimize()
    new_fba = FBA(new_model).build().optimize()

    try:
        assert round_(fba.objective_value) == round_(new_fba.objective_value)
    except:
        print(colored("FBA objective values differ", "red"))
        print(colored(f"  new fba: {new_fba.objective_value}", "blue"))
        print(colored(f"  orig fba: {fba.objective_value}", "blue"))
        set_output_value(FAILED)


def do_fva_asserts(original_model, new_model, fraction=1.0, reactions=False):

    if reactions:
        reaction_list = list(original_model.reactions.keys())[0:40]
    else:
        reaction_list = None

    orig_fva = fva(original_model, fraction=fraction, reactions=reaction_list)
    new_fva = fva(new_model, fraction=fraction, reactions=reaction_list)
    
    try:
        assert compare_fva_dataframes(new_fva, orig_fva)
    except:
        set_output_value(FAILED)


def do_srfba_asserts(original_model, new_model):
    srfba = SRFBA(original_model).build().optimize()
    new_srfba = SRFBA(new_model).build().optimize()

    try:
        assert round_(srfba.objective_value) == round_(new_srfba.objective_value)
    except:
        print(colored("SRFBA objective values differ", "red"))
        print(colored(f"  new srfba: {new_srfba.objective_value}", "blue"))
        print(colored(f"  orig srfba: {srfba.objective_value}", "blue"))
        set_output_value(FAILED)


def do_rfba_asserts(original_model, new_model):
    rfba = RFBA(original_model).build().optimize()
    new_rfba = RFBA(new_model).build().optimize()

    try:
        assert round_(rfba.objective_value) == round_(new_rfba.objective_value)
    except:
        print(colored("RFBA objective values differ", "red"))
        print(colored(f"  new rfba: {new_rfba.objective_value}", "blue"))
        print(colored(f"  orig rfba: {rfba.objective_value}", "blue"))
        set_output_value(FAILED)


if __name__ == "__main__":


    if e_coli:
        # Tests for E Coli models
        e_coli_gem_file = "/home/sebastiao_zoio/Documents/test_mewpy/mewpy/examples/models/germ/e_coli_core.xml"

        e_coli_trn_model = "/home/sebastiao_zoio/Documents/test_mewpy/mewpy/examples/models/germ/e_coli_core_trn.csv"
     
        #e_coli_original_germ_model = get_orginal_gem_model(e_coli_gem_file)
        #e_coli_original_germ_model = get_orginal_gem_model(e_coli_gem_file, e_coli_trn_model)
        e_coli_new_germ_model = get_new_gem_model(e_coli_gem_file, e_coli_trn_model)

        #do_fba_asserts(e_coli_original_germ_model, e_coli_new_germ_model)

        #do_fva_asserts(e_coli_original_germ_model, e_coli_new_germ_model, reactions=True)

        #do_srfba_asserts(e_coli_original_germ_model, e_coli_new_germ_model)

        #do_rfba_asserts(e_coli_original_germ_model, e_coli_new_germ_model)


    if yeast:
        # Tests for Yeast models
        yeast_gem_file = "/home/sebastiao_zoio/Documents/test_mewpy/mewpy/examples/models/germ/yeast-GEM.xml"

        yeast_original_germ_model = get_orginal_gem_model(yeast_gem_file)
        yeast_new_germ_model = get_new_gem_model(yeast_gem_file)

        do_fba_asserts(yeast_original_germ_model, yeast_new_germ_model)

        do_srfba_asserts(yeast_original_germ_model, yeast_new_germ_model)

        do_rfba_asserts(yeast_original_germ_model, yeast_new_germ_model)

    if output_value:
        print(colored("Success", "green"))