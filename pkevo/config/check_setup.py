from sys import version_info
import os
import inspect
# Note: rdkit is not needed here, but importing it now prevents possible Boost.Python errors later on
from rdkit import DataStructs
from rdkit import Chem
from rdkit.Chem import AllChem

from config.pkevo_config import SearchConfig as search_config

# Define the class_names list as a global variable
fitness_function_classnames = []

# the first of 2 steps
def initial_setup():
    """
    Performs the following checks and raises exceptions or adjusts the setup:
        - check Python version > 3.8 (throws error if not fulfilled)
        - check if logs folder exists (initializes FileStorage Singleton and creates logs folder if necessary)
        - check which classes are available in grammatical_evolution.ge_utils.fitness_functions.py
    """

    # check requirede minimal Python version (3.8 required for f-strings new '{var=}' feature) 
    if version_info.major < 3 or (version_info.minor < 8 and version_info.major == 3):
        print("\nError: Unsupported Python version! Please use Python 3.0 or higher")
        exit(1)

    # Initialize FileStorage class and check if logs folder exists
    from config.file_storage_singleton import file_storage

    if not os.path.exists(file_storage.logs_dir):
        file_storage.create_logs_folder()

    # Check which classes are defined in "fitness_functions.py"
    try:
        from grammatical_evolution.ge_utils import fitness_functions

        # Get all classes defined in "fitness_functions.py"
        classes = inspect.getmembers(fitness_functions.FitnessFunctions, inspect.isclass)

        # Get the names of all Fitness classes (ignore template and __class__) and store them in the global variable
        global fitness_function_classnames
        fitness_function_classnames = [cls[0] for cls in classes if cls[0] != "__class__" and cls[0] != "FunctionTemplate"]
    except Exception as e:
        print(f"An error occurred while checking 'fitness_functions.py': {str(e)}")
        exit(1)


def _check_fitness_class_attributes(fitness_class):
    # List of required attributes and methods
    required_attributes = ["coevolution", 
                           "default_fitness", 
                           "optimization_type", 
                           "is_chemistry_necessary", 
                           "target", 
                           "objective_function"]

    # Check if the class has the required attributes
    for attr in required_attributes:
        if not hasattr(fitness_class, attr):
            print(f"AttributeError: The Fitness function {fitness_class.__name__} is missing the following: {attr}")
            exit(1)


# second and final step
def adjust_setup(input_args):
    """
    Performs the following checks:
        1) check if simulation_results folder is required (creates folder if necessary)
        2) check if selected fitness function is valid
        3) check if all requirements for selected fitness function are fulfilled
        4) check if user provided optional arguments, which would require the configuration parameters to be adjusted
    If case of problems, setup is either adjusted or the app is exited with an error.
    """
    from grammatical_evolution.ge_utils.fitness_functions import FitnessFunctions
    import config.custom_logger as custom_logger
    logger = custom_logger.get_logger(__name__, log_file='config_logs.log')
    # Note: from now on forwards we can use the custom logger in the application


    # Check 1: Unless user explicitly stated not to create any simulation result files, create a new simulation folder
    if not input_args.no_files:
        from config.file_storage_singleton import file_storage
        file_storage.create_simulation_folder()


    # Check 2: see if selected fitness function is available
    fitness_function_name = input_args.fitness_function

    if input_args.fitness_function not in fitness_function_classnames:
        # Check if user provided only starting letters of desired fitness function
        found_alias = False
        for ff_classname in fitness_function_classnames:
            if ff_classname.lower().startswith(input_args.fitness_function.lower(), ):
                input_args.fitness_function = ff_classname
                fitness_function_name = ff_classname
                found_alias = True
                break
        if not found_alias:    
            print(f"ValueError: Unknown fitness function '{input_args.fitness_function}'. Available options are {fitness_function_classnames}")
            exit(1)

    # Try loading the fitness function class and check if class has all requried attributes/methods
    try:
        selected_fitness_class = getattr(FitnessFunctions, fitness_function_name)
        _check_fitness_class_attributes(selected_fitness_class)
    except AttributeError as e:
        print(e)
        exit(1)
    except Exception as e:
        print(f"An error occurred while loading the fitness function {fitness_function_name}: {str(e)}")
        exit(1)


    # Check 3: see if requirements for fitness function are fulfilled

    # For the pharmacophore fitness function, we need both MØD and LigandScout
    if fitness_function_name == 'PharmacophoreScore':
        # Try to import LigandScout
        from external_tools.ligand_scout_singleton import ligand_scout_utils
        ligang_scout_version = ligand_scout_utils.get_ligand_scout_version()
        if ligang_scout_version == 'unknown':
            logger.error("Cannot use pharmacophore fitness for evaluation, because LigandScout is not available. Using 'FingerprintSimilarity' instead.")
        # Try to import MØD
        try:
            from external_tools.mod_importer import mod
        except ImportError:
            logger.error("Cannot use pharmacophore fitness for evaluation, because MOD is not available. Using 'StringDistance' instead.")

    # For the fingerprint similarity of molecules we need MØD
    elif fitness_function_name == 'FingerprintSimilarity':
        try:
            from external_tools.mod_importer import mod
        except ImportError:
            logger.error("Cannot use molecular fingerprint similarity for evaluation, because MOD is not available. Using 'StringDistance' instead.")


    # Check 4: Look for optional parameters passed by the user, which would require our pkevo_config parameters to be updated
    # Note: As also mentioned in the wiki, arguments passed by the user via CLI arguments take precedence over what is defined in pkevo_config.py
    if input_args.threshold:
        search_config.FITNESS_THRESHOLD = input_args.threshold
    if input_args.break_after:
        search_config.BREAK_AFTER = input_args.break_after

    return input_args
