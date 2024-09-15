import inspect

from grammatical_evolution.ge_utils.fitness_functions import FitnessFunctions
from config.check_setup import _check_fitness_class_attributes

def test_class_definitions():

    # Load all defined fitness functions from "fitness_functions.py"
    fitness_function_classnames = []
    try:
        from grammatical_evolution.ge_utils import fitness_functions

        # Get all classes defined in "fitness_functions.py"
        classes = inspect.getmembers(fitness_functions.FitnessFunctions, inspect.isclass)

        # Get the names of all Fitness classes (ignore template and __class__) and store them in the global variable
        fitness_function_classnames = [cls[0] for cls in classes if cls[0] != "__class__" and cls[0] != "FunctionTemplate"]

        # Iterate over available fitness classes and check if they have all required attributes defined
        for fitness_function_name in fitness_function_classnames:
            selected_fitness_class = getattr(FitnessFunctions, fitness_function_name)
            # Note: I know, don't call private functions of other modules in your production code! But IMO, for a test case it's acceptable
            _check_fitness_class_attributes(selected_fitness_class)
    except Exception as e:
        exit(e)
