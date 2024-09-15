import argparse
import os

from config.pkevo_config import EvolutionConfig as evo_config
from config.check_setup import fitness_function_classnames

def parse_args():
    config_file_path = os.path.join(os.path.dirname(__file__), "pkevo_config.py")
    avail_fit_functions = fitness_function_classnames
    # Set default fitness function to Pharmacophore Score if available, otherwise pick first one from list
    default_fit_function = 'PharmacophoreScore' if 'PharmacophoreScore' in avail_fit_functions else avail_fit_functions[0]

    parser = argparse.ArgumentParser(
        usage=argparse.SUPPRESS,
        description=f"""<========= Welcome to the PKevo Help page =========>
            `PK` stands for `P`oly`K`etide (and definitely not for the initials of the developer). The app lets you run a GE 
            (= Grammatical Evolution) simulation for Polyketides. You can configure parameters in the '{config_file_path}' file.""",
        epilog="To run a GE simulation execute in your CLI the command 'python main.py' (appended with above described arguments if desired)")
    
    parser.add_argument("-v", "--verbose", 
                        action="store_true", 
                        help="Use this flag to get additional output printed to the CLI")
    
    parser.add_argument("-si", "--use-starting-inds", 
                        action="store_true", 
                        help="This flag enables initialization of the GE with the individuals defined in /config/starting_individuals.json")
    
    parser.add_argument("-nf", "--no-files", 
                        action="store_true", 
                        help="Use this flag to disable file generation for storing results and other data")
    
    parser.add_argument("-f", "--fitness-function",  
                        default=default_fit_function,
                        help=f"Choose the fitness function to use: {avail_fit_functions}. Defaults to: '{default_fit_function}'")

    parser.add_argument("-t", "--threshold",
                        type=float,  
                        help="If set (as a float value), the GE algorithm will stop once the defined fitness threshold is surpassed")
    
    parser.add_argument("-b", "--break-after",
                        type=int,
                        help="If set (as an int value), the GE algorithm will stop if the overall best fitness has not improved for x GE iterations")
    
    return parser.parse_args()