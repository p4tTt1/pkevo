import os

# this initial setup needs to happen before anything else!
from config.check_setup import initial_setup, adjust_setup
initial_setup()

import config.argument_parser as parser
import config.custom_logger as custom_logger
from config.file_storage_singleton import file_storage
from config.pkevo_config import EvolutionConfig as evo_config

from grammatical_evolution.grammatical_evolution import GrammaticalEvolution
#from grammatical_evolution.ge_utils.fitness_functions import FitnessFunctions

#import logging
#logging.disable() # Note: Fast and easy way to disable all logging
logger = custom_logger.get_logger(__name__, log_file='ge_logs.log')


def main():
    logger.info("Started a new simulation.")
    # Parse arguments received from CLI
    args = parser.parse_args()

    # Check the current app setup against the provided user arguments
    adjust_setup(args)

    # Get the current directory where this script is located
    current_dir = os.path.dirname(os.path.abspath(__file__))

    # Construct the path to the grammar file
    grammar_file_path = os.path.join(current_dir, 'grammatical_evolution', 'grammars', evo_config.GRAMMAR_FILE_NAME)

    # Create the GrammaticalEvolution object
    ge = GrammaticalEvolution(grammar_file_path, 
                              fitness_function_name=args.fitness_function,
                              verbose=args.verbose, 
                              store_results= not args.no_files, 
                              use_starting_inds=args.use_starting_inds)

    # Run the Grammatical Evolution algorithm
    best_ever = ge.run(verbose=args.verbose)

    # Print the best individual found
    print("\n=======================")
    print("Best individual:")
    print(best_ever.__str__(verbosity='vv'))
    print("=======================\n")


# Ensure the main code is only executed when ran as a script and not when imported by someone
if __name__ == "__main__":
    main()
