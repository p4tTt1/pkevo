""" 
Based on Robert Haas' "panakeias_garden" package (/panakeias_garden/pre_release_demos/metaheuristics.py) 
and PonyGE version 1 (https://github.com/jmmcd/ponyge).

The original PonyGE was chosen over the newer PonyGE2 as a starting-point, since it is written in a single file. 
This makes it easier to adapt it to the specific needs of PKevo, instead of PKevo having to adapt to the design 
choices made by PonyGE2.
"""

import copy
import random
import time
from typing import List
from multiprocessing import Pool

from config.pkevo_config import EvolutionConfig as evo_config
from config.pkevo_config import SearchConfig as search_config
from config.pkevo_config import GeneralConfig
import config.custom_logger as custom_logger
from config.file_storage_singleton import file_storage

from generation_management.generation_snapshots import GenerationSnapshots

from .ge_utils.fitness_functions import FitnessFunctions
from .ge_utils import selection_strategies, replacement_strategies, mutations, crossovers
from .ge_utils import stats_calculator as stats_calculator
from .ge_utils.evolution_plotter import EvolutionPlotter

from .model._grammar import _Grammar
from .model._individual import _Individual
from .model._molecule import Molecule

logger = custom_logger.get_logger(__name__, "ge_logs.log")


class GrammaticalEvolution():
    """
    Represents the Grammatical Evolution algorithm for evolving individuals according to a context-free grammar.
    Main attributes:
        grammar (_Grammar): An instance of the _Grammar class representing the grammar used for the Grammatical evolution
        fitness_function (function): The fitness function used to evaluate individuals.
        chemistry (Chemistry or None): An instance of the Chemistry class, used for specific fitness functions.
        file_storage (FileStorageSingleton or None): An instance of FileStorageSingleton for storing results persistently as files.
        starting_inds (list or None): optional list containing genomes and phenotypes for individuals to be used for "Gen 0".
    Main methods:
        run(verbose=False):
            Runs the Grammatical Evolution algorithm, including initialization, evolution loop, and result storage.
        evo_loop(max_generations, individuals, verbose=False):
            Takes care of looping over all evolutionary generations. The actual evolution is done in the evo_step() function however!.
        evo_step(individuals, best_ever, idx):
            Performs one step of the evolutionary process, including selection, crossover, mutation, and replacement.
            Returns the updated population and the best individual found during the step.
    """
    

    def __init__(self, grammar_file_path, fitness_function_name, verbose=False, store_results=True, use_starting_inds=False):
        """
        Constructor for the GrammaticalEvolution class.
        Args:
            grammar_file_path (str): The file path to the grammar used in evolution.
            fitness_instance (FitnessFunctions.FunctionTemplate): The selected Fitness function class
            verbose (bool, optional): Whether to enable verbose mode with additional output (default is False).
            store_results (bool, optional): Whether to enable result storage (default is True).
            use_starting_inds (bool, optional): if true, GE will try to use defined individuals for initial population (default is False).
        """
        self.grammar = _Grammar(grammar_file_path=grammar_file_path, verbose=verbose, store_results=store_results)
        self.gen_snapshots = GenerationSnapshots()

        self.remembered_fitness_scores = dict() # Initialize an empty dictionary for remembering fitness scores
        self.unique_genomes = dict() # dictionary with `genome_tuple` (just the genome list turned into a tuple) as key and _Individual obj as value
        self.counter_total_ind = 0 # for debugging only
        self.counter_eval_fit = 0 # for debugging only
        self.counter_released_mols = 0 # for debugging only

        # Dynamically load the desired Fitness function class from fitness_functions.py (no checks needed, already done in `check_setup` script)
        selected_fitness_class = getattr(FitnessFunctions, fitness_function_name)
        # Create an instance of the selected class
        fitness_function_instance: FitnessFunctions.FunctionTemplate = selected_fitness_class()
        self.fitness_function = fitness_function_instance

        # Select operators and strategies (dynamically loading the methods specified in the config file)
        self.replacement_strategy = getattr(replacement_strategies, evo_config.REPLACEMENT_STRATEGY)
        self.selection_strategy = getattr(selection_strategies, evo_config.SELECTION_STRATEGY)
        self.crossover_operator = getattr(crossovers, evo_config.CROSSOVER_OPERATOR)

        if store_results:
            self.file_storage = file_storage
            self.file_storage.save_parameters(evo_config, search_config) # Creates a file storing all used parameters for this GE run
        else:
            self.file_storage = None

        if use_starting_inds:    
            self.starting_inds = file_storage.prepare_individuals()
        else:
            self.starting_inds = None

        # Setup the chemistry space (available molecules, reactions and substrategies).
        # This is only to be done once during the "initialization" of the whole process.
        # Later on only the specific strategy for each given PKS Domain string will be 
        # generated on-the-fly and used for the derivation graph construction.
        if fitness_function_instance.is_chemistry_necessary:
            from chemistry.chemistry import Chemistry 
            chemistry = Chemistry()
            self.chemistry = chemistry
        else:
            # There is no need to get chemistry involved
            self.chemistry = None

        # Create an instance of the EvolutionPlotter class
        self.plotter = EvolutionPlotter(fitness_function_instance.optimization_type, evo_config.GENERATIONS)
    

    def run(self, verbose=False)->_Individual:
        """
        Initializes the population, selects the genetic operators (crossover and mutation), and runs the loop over all remaining generations.
        Args:
            verbose (bool, optional): Whether to enable verbose mode with additional output (default is False).
        Returns:
            Individual: The best Individual evolved.
        """
        # Create Individuals
        individuals = self.initialize_population(evo_config.POPULATION_SIZE, self.starting_inds)
        # Note: at this point, an individual has a genotype and phenotype, but no fitness yet
        
        # Let's go into the evolutionary loop...
        print("=======================")
        print("Starting evolutionary loop...")
        print("=======================")
    
        best_ever, duration = self.evo_loop(evo_config.GENERATIONS, individuals, verbose)

        print("======DEBUG========")
        print(f"{len(self.unique_genomes)=}")
        print(f"{self.counter_total_ind=}")
        print(f"{self.counter_eval_fit=}")
        if self.chemistry:
            print(f"{self.counter_released_mols=}")
        print("======DEBUG========")

        print(f"\nSimulation finished after {duration:.3f} seconds.")
        return best_ever


    def evo_loop(self, max_generations, individuals: List[_Individual], verbose=False):
        """
        The main evolutionary loop that performs selection, crossover, mutation, and fitness evaluation. During that it
        also takes care of generating statistics for later evaluation of the simulation.
        Args:
            max_generations (int): The maximum number of generations to run the evolutionary loop.
            individuals (List[_Individual]): A list of individuals representing the current population.
            verbose (bool, optional): Whether to enable verbose mode with additional output (default is False).
        Returns:
            Individual: The best individual evolved during the GE algorithm.
            duration: The duration of the evolutionary loop in seconds
        """
        start_time_loop = time.perf_counter()
        all_generation_info = []
        if self.starting_inds is not None:
            starting_info = f"Included {len(self.starting_inds)} starting individuals in initial population.\n"
        else:
            starting_info = "No specific starting individuals were used.\n"
        all_generation_info.append(starting_info)

        # Loop through all generations
        best_ever = None
        current_best_fitness = self.fitness_function.default_fitness
        runs_without_improvement = 0
        generation_count = 0           

        for generation_idx in range(0, max_generations):
            start_time_step = time.perf_counter()
            individuals, best_ever = self.evo_step(individuals, best_ever, generation_idx)
            end_time_step = time.perf_counter()
            elapsed_time_step = end_time_step - start_time_step
            individuals[0].build_pks_architecture() # build the architecture of best ind in this generation (interesting for result analysis)
            fitness_stats, used_codons_stats, phenotype_stats = stats_calculator.calculate_statistics(individuals)
            generation_count = generation_idx + 1

            # Update current best fitness if it is better than previous one (needed for break condition)
            if best_ever.fitness != current_best_fitness:
                # We found a new, better fitness -> update the current best fitness value and reset `runs_without_improvement`
                current_best_fitness = best_ever.fitness
                runs_without_improvement = 0
            else:
                # We did not improve in the current iterations -> increment `runs_without_improvement`
                runs_without_improvement += 1
            
            # Generate string containing information about the current generation
            generation_info = (
                f"\nGeneration #{generation_idx} (Elapsed time: {elapsed_time_step:.3f}s): Fitness min:{fitness_stats[0]:.2f} avg:{fitness_stats[2]:.2f} max:{fitness_stats[1]:.2f}"
                f"\nPhenotypes unique: {phenotype_stats[0]}, duplicated: {phenotype_stats[1]}, invalid: {phenotype_stats[2]}"
                f"\nBest {individuals[0]}"
            )

            # Append to all_generation_info if result storage is enabled
            if self.file_storage:
                all_generation_info.append(generation_info)        

            # Print generation info if verbose is enabled
            if verbose:
                print(generation_info)
            else:
                # show a progress bar in CLI instead (otherwise with no feedback user might think the app got stuck)
                self._visualize_progress(generation_idx+1, elapsed_time_step)        

            # If plotting is enabled, add the data for current generation to plot data
            if GeneralConfig.GENERATE_PLOTS:
                self.plotter.add_data(fitness_stats, used_codons_stats, phenotype_stats, generation_idx)
            
            # Check break conditions to see whether or not we should continue with the evolutionary loop
            need_a_break = self._is_break_condition_reached(best_ever.fitness, runs_without_improvement)
            if need_a_break:
                # At least one of the break conditions has been reached, so we should quit the loop now
                break
        # End of evolutionary loop => stop the timer    
        end_time_loop = time.perf_counter()
        elapsed_time_loop = end_time_loop - start_time_loop
        
        # Create some plots for easier interpretation of the results, unless disabled by user choice
        if GeneralConfig.GENERATE_PLOTS:
            self.plotter.plot_evolution(generation_count, self.file_storage)

        # Optionally create a file showing the evolution over the generations
        if self.file_storage:
            # Store best individual
            file_storage.save_best_ind(best_ever)
            # Create image file for best molecule, if molecules exist
            if len(best_ever.molecules) > 0:
                file_storage.save_molecule(best_ever.molecules[0].smiles_string)
            # Store data about all generations
            all_generation_info.insert(0, (f"Simulation finished after {elapsed_time_loop:.3f} seconds."))
            simulation_info = '\n'.join(all_generation_info)
            self.file_storage.save_text_file("generations.txt", simulation_info)
            # Store the final generation as a snapshot
            self.gen_snapshots.add_generation(iteration = generation_count, generation = individuals)
            self.gen_snapshots.export_snapshots()

        return best_ever, elapsed_time_loop


    def evo_step(self, curr_pop, best_ever, loop_count):
        """
        Perform one step of the evolutionary process, including selection, crossover, mutation, and fitness evaluation.
        Args:
            curr_pop (List[_Individual]): The current population of individuals.
            best_ever (_Individual): The best individual evolved up to the current generation.
            loop_count (int): The current generation number.
        Returns:
            tuple: A tuple containing the new population of individuals and the best individual evolved in this step.
        """
        new_pop = []

        if loop_count == 0:
            # For Init Generation no selection, crossover, mutations and replacement, however we need to sort the population and select best_ever  manually
            self.perform_evaluation(curr_pop)
            curr_pop.sort(reverse=True)
            new_pop = curr_pop
            best_ever = new_pop[0]
            # Store the initial generation as a snapshot
            if self.file_storage:
                self.gen_snapshots.add_generation(iteration = 0, generation = new_pop)
        else:
            # Select parents
            parents = self.selection_strategy(curr_pop)
            # Crossover parents and add offspring to the new population
            while len(new_pop) < evo_config.GENERATION_SIZE:
                # 2 new individuals are created by selecting 2 parents randomly and performing the crossover operation
                new_pop.extend(self.perform_crossover(*random.sample(parents, 2), evo_config.WITHIN_CODING_REGION, evo_config.VARIABLE_LENGTH))
            # Mutate the new population and give them their phenotypes
            new_pop = [self.get_phenotype(self.perform_mutation(new_ind, evo_config.WITHIN_CODING_REGION)) for new_ind in new_pop]
        
            # Evaluate the fitness of the new population
            self.perform_evaluation(new_pop)
         
            # Replace the current population (already sorted) with the new population (not yet sorted!) depending on selected replacement strategy
            new_pop = self.replacement_strategy(new_pop, curr_pop) # replacing also sorts the population accordingly
            best_ever = max(best_ever, new_pop[0]) # Q: isn't best_ever garantueed to be part of individuals? -> Not necessarily, depends on replacement strategy

        return new_pop, best_ever


    def _is_break_condition_reached(self, curr_best_fitness, runs_without_improvement):
        """
        This function determines whether to exit the evolutionary loop early or not based on defined thresholds.
        Arguments:
            curr_best_fitness (float): The overall best fitness so far achieved by any of the generated individuals
            runs_without_improvement (int): The number of evolutionary steps since the best overall fitness has been improved
        Return:
            boolean: True if a break condition was satisfied, False otherwise
        """
        if search_config.FITNESS_THRESHOLD is not None:
            if (self.fitness_function.optimization_type == 'min' and curr_best_fitness <= search_config.FITNESS_THRESHOLD) or \
               (self.fitness_function.optimization_type == 'max' and curr_best_fitness >= search_config.FITNESS_THRESHOLD):
                logger.info(f"Exiting evolutionary loop because we reached threshold with this individual's fitness: {curr_best_fitness}")
                print(f"Exiting evolutionary loop because we reached threshold with this individual's fitness: {curr_best_fitness}")
                return True
        if search_config.BREAK_AFTER:
            if runs_without_improvement >= search_config.BREAK_AFTER:
                logger.info(f"Exiting evolutionary loop early, since we haven't improved in {runs_without_improvement} generations")
                print(f"Exiting evolutionary loop early, since we haven't improved in {runs_without_improvement} generations")
                return True

        return False 


    def _visualize_progress(self, count, duration):
        """
        Helper function to show progress in CLI (useful when 'verbose' is set to false, yet we want to indicate to the user
        that the algorithm is still working and not stuck somewhere.
        Args:
            count (int): The number of the current generation.
            duration (float): Time (in seconds) it took to complete the evolutionary step for this generation
    
        From https://stackoverflow.com/questions/3173320/text-progress-bar-in-terminal-with-block-characters#answer-27871113
        """
        percents = round(100.0 * count / float(evo_config.GENERATIONS), 1)
        print(f"Simulated generation {count} out of {evo_config.GENERATIONS} in {duration:.3f} seconds. ({percents}% completed)")


    def initialize_population(self, size=10, starting_inds=None):
        """
        Initialize the population of individuals.
        Args:
            size (int, optional): The size of the population (default is 10).
            starting_inds (list, optional): list containing genomes and phenotypes of individuals (default is None).
        Returns:
            List[_Individual]: A list of initialized individuals representing the population.
        """
        population = []

        if starting_inds is not None:
            for starting_ind in starting_inds:
                # In case the population has already reached the set population size, we must aboard the process
                if len(population) == size:
                    return population
                # Check that the given Genome is able to produce the desired phenotype with the current grammar
                start_phenotype = self.grammar.generate_phenotype(starting_ind[0])
                if start_phenotype[0] == starting_ind[1]:
                    # We got a match! Desired starting individual can be added to initial population
                    logger.info("Using the desired starting individual with phenotype: " + start_phenotype[0])
                    population.append(self.get_phenotype(_Individual(self.fitness_function.optimization_type, self.fitness_function.default_fitness, genome=starting_ind[0])))
                else:
                    logger.warning(f"Ignoring the desired starting individual ({start_phenotype[0]}) since it could not be generated with the current grammar.")
            # If we weren't able to generate any valid starting individuals, unset the starting_inds parameter
            if len(population) == 0:
                self.starting_inds = None
        # Fill the rest of the population with randomly generated individuals and give them their phenotypes
        population.extend([self.get_phenotype(_Individual(self.fitness_function.optimization_type, self.fitness_function.default_fitness)) for _ in range(size - len(population))])

        return population


    def create_random_ind(self)-> _Individual:
        """
        Create a random individual (for testing and debugging purposes).
        Returns:
            _Individual: A randomly generated individual.
        """
        rand_ind = self.initialize_population(1)[0] # init_pop returns an array, so let's make sure we return only the single element from the list
        rand_ind = self.get_phenotype(rand_ind)
        return rand_ind


    def perform_mutation(self, individual: _Individual, within_used=True):
        """
        Perform different types of mutations on an individual's genotype with a predefined probability.
        Args:
            individual (_Individual): The individual to mutate.
            within_used (bool, optional): Whether to perform mutations only within used codons region (default is True).
        Returns:
            _Individual: The mutated individual.
        """
        # Perform single-point mutations of the genome
        mutations.int_flip_mutation(individual=individual, within_used=within_used)

        # Perform deletion mutation with predefined probability
        if random.random() < evo_config.DELETION_PROBABILITY:
            mutations.deletion_mutation(individual=individual, within_used=within_used)

        # Perform insertion mutation with predefined probability
        if random.random() < evo_config.INSERTION_PROBABILITY:
            mutations.insertion_mutation(individual=individual, within_used=within_used)

        return individual


    def perform_crossover(self, p_0: _Individual, p_1: _Individual, within_used=True, variable_len=True):
        """
        Given two individuals (parents), create two children using the predefined crossover operation and return them.
        Args:
            p_0 (_Individual): The first parent individual.
            p_1 (_Individual): The second parent individual.
            within_used (bool, optional): Whether to perform crossover only within used codons region (default is True).
            variable_len (bool, optional): If True, resulting offspring genomes might vary in lenght (Defaults to True.)
        Returns:
            Tuple[_Individual, _Individual]: A tuple containing two child individuals resulting from crossover.
        """
        return self.crossover_operator(p_0, p_1, self.fitness_function.default_fitness, within_used, variable_len)


    def get_phenotype(self, ind: _Individual):
        """
        This function gives an individual its phenotype by either performing a call to `Grammar.generate_phenotype()` 
        for previously unknown individuals (those with unique genomes), which sets the phenotype, used_codons and the 
        assembly_line for an individual. Already known individuals get copied.
        Args:
            ind (_Individual): The individual for which the phenotype should be generated
        Returns:
            _Individual: The updated individual
        """
        # Let's check if the current individual already exists in our dict() cache, so that we don't run calculations repeatedly
        genome_tuple = tuple(ind.genome)  # Convert the list to a tuple (so that we can use it as an index)
        if genome_tuple not in self.unique_genomes:
            ind.phenotype, ind.used_codons, ind.assembly_line = self.grammar.generate_phenotype(ind.genome) # Decoding step: Genotype -> Phenotype
        else:
            # load the individual from the cached dictionary
            ind = copy.copy(self.unique_genomes[genome_tuple])

        return ind


    def get_molecules(self, ind: _Individual):
        """
        Creates molecules using the `Chemistry` package for an individual if it has a PKS assembly line and it doesn't
        already contain molecules from past calls.
        Args:
            ind (_Individual): The Individual for which molecules should be generated
        Returns:
            List[Molecule]: A list of Molecule objects, which were generated for this Individual 
        """
        # In case the individual was retrieved from cache, it might already have molecules, which we can return directly
        if len(ind.molecules) > 0:
            return ind.molecules

        molecules = []
        # Check if Individual is valid, otherwise it cannot generate molecules
        if ind.assembly_line:
            molecules = [Molecule(smiles_string, 0.0) for smiles_string in self.chemistry.get_smiles_for_ind(ind.assembly_line)]
 
        return molecules


    def perform_evaluation(self, individuals: List[_Individual]):
        """
        This function takes care of evaluating the fitness of a list of individuals by: 
            1) If necessary, triggering the creation of molecules for the Individuals
            2) Calling the user specified fitness function for either the whole population or each individual separately
            3) If molecules are at play, the list of molecules is sorted by their fitness
            4) Storing fitness results in the self.unique_genomes dictionary for future lookup
        Args:
            (List[_Individual]): A list of individuals to evaluate.
        """
        # Get chemistry involved if necessary and generate molecules for all valid individuals.
        if self.chemistry:
            # Use Multiprocessing to speed things up
            if GeneralConfig.MULTIPROCESSING_CORES > 1:
                # Set an appropriate chunksize (which is the number of tasks per process), commonly known as batch size
                # Doing this will initialize the Chemistry class for each process, but that takes only a fraction of a second
                batch_size = int(len(individuals) / GeneralConfig.MULTIPROCESSING_CORES) + 1
                with Pool(processes=GeneralConfig.MULTIPROCESSING_CORES) as pool:
                    molecules = pool.map(func=self.get_molecules, iterable=individuals, chunksize=batch_size)
                    
                # Update the individuals list with the newly generated molecules    
                for ind, molecule_list in zip(individuals, molecules):
                    self.counter_released_mols += len(molecule_list)
                    ind.molecules = molecule_list
            else:
                # No multiprocessing
                for ind in individuals:
                    ind.molecules = self.get_molecules(ind)
                    self.counter_released_mols += len(ind.molecules)

        # Evaluate fitness of the population as a whole
        if self.fitness_function.coevolution:
            self._evaluate_generation(individuals)
        else:
            # Iterate over population and evaluate each individual separately
            for ind in individuals:
                self._evaluate_individual(ind)


    def _evaluate_generation(self, individuals: List[_Individual]):
        """
        This evaluates the whole generation at once, allowing for co-evolutionary factors to be 
        taken into account. 
        
        Note: Since co-evolution might be involved here, it is not possible to load previous fitness scores
        for 'identical' individuals from previous evaluations, since due to co-evolution the same phenotype might
        be scored differently this time.

        Args:
            (List[_Individual]): A list of individuals to evaluate.
        """
        # Evaluate whole generation instead of each individual separately
        self.fitness_function.objective_function(individuals, self.remembered_fitness_scores)
        
        # Nonetheless, we need to iterate over all individuals once again to perform molecule sorting and storing Inds to our cache
        for ind in individuals:
            if len(ind.molecules) > 1:
                ind.sort_molecules()
            genome_tuple = tuple(ind.genome)  # Convert the list to a tuple (so that we can use it as an index)
            if genome_tuple not in self.unique_genomes: 
                self.unique_genomes[genome_tuple] = ind


    def _evaluate_individual(self, ind: _Individual):
        """
        Counter part to the `_evaluate_generation` function. Here only a single individual is evaluated.
        
        Args:
            _Individual: An individual to evaluate.
        """
        self.counter_total_ind += 1

        # Let's check if the current individual already exists in our dict() cache, so that we don't run calculations repeatedly
        genome_tuple = tuple(ind.genome)  # Convert the list to a tuple (so that we can use it as an index)
        if genome_tuple not in self.unique_genomes:
            # We can only evaluate individuals who have a phenotype, otherwise we leave them with their default fitness
            if ind.phenotype is not None:
                self.counter_eval_fit += 1
                ind.fitness = self.fitness_function.objective_function(ind, self.remembered_fitness_scores)                   
                # If Individual has multiple molecules, make sure this list is sorted properly
                ind.sort_molecules()
            # Finally add the new and evaluated individual to our "Individuals database"
            self.unique_genomes[genome_tuple] = ind
        else:
            # An Individual with the same genome was evaluated once before -> let's copy its fitness and potential molecules
            ind.fitness = self.unique_genomes[genome_tuple].fitness
            ind.molecules = copy.copy(self.unique_genomes[genome_tuple].molecules)
