import os
import pytest
import copy
import json

from grammatical_evolution.grammatical_evolution import GrammaticalEvolution
from grammatical_evolution.model._individual import _Individual
from grammatical_evolution.ge_utils.fitness_functions import FitnessFunctions

from config.pkevo_config import EvolutionConfig as evo_config
from config.pkevo_config import SearchConfig as search_config


@pytest.fixture
def ge_instance():
    # Get the current directory where this script is located
    current_dir = os.path.dirname(os.path.abspath(__file__))

    # Construct the path to the grammar file
    grammar_file_path = os.path.join(current_dir, '..', 'grammatical_evolution', 'grammars', 'pytest_pks.bnf')

    # Create the GrammaticalEvolution object
    ge = GrammaticalEvolution(grammar_file_path, 
                              fitness_function_name="StringDistance",
                              verbose=False, 
                              store_results=False, 
                              use_starting_inds=False)
    return ge


def load_json_data():
    current_dir = os.path.dirname(os.path.abspath(__file__))
    json_file_path = os.path.join(current_dir, 'starting_inds_test.json')

    with open(json_file_path, 'r') as file:
        data = json.load(file)

    individuals = []
    for entry in data:
        genome = entry.get("Genome", [])
        phenotype = entry.get("Phenotype", "")
        individual = [genome, phenotype]
        individuals.append(individual)

    return individuals


# Basic test to see if an instance of the GrammaticalEvolution can be initialized properly
def test_ge_init(ge_instance: GrammaticalEvolution):
    assert ge_instance is not None
    assert ge_instance.grammar is not None
    assert ge_instance.fitness_function is not None
    assert ge_instance.chemistry is None
    assert ge_instance.starting_inds is None
    assert ge_instance.file_storage is None
    assert ge_instance.counter_total_ind == 0


## Test for population initialization (with no starting ind)
def test_initialize_population(ge_instance: GrammaticalEvolution):
    # Initialize a new population with all random individuals
    pop = ge_instance.initialize_population(size=10, starting_inds=None)
  
    assert isinstance(pop, list)
    assert len(pop) == 10
    assert isinstance(pop[0], _Individual)
    
    # Let's check the properties of an Individual
    ind_0: _Individual = pop[0]
    assert len(ind_0.genome) == 100
    assert isinstance(ind_0.phenotype, str) or isinstance(ind_0.phenotype, None)
    assert ind_0.used_codons >= 0
    # Architecture is defined at a later stage of the evolutionary loop
    assert ind_0.architecture is None
    # Fitness has not yet been evaluated and should therefore remain at its default value
    assert ind_0.fitness == FitnessFunctions.StringDistance.default_fitness


## Test for population initialization (with valid starting ind)
def test_init_pop_valid_starter_ind(ge_instance: GrammaticalEvolution):
    # initialize population with a valid individuals (valid means that provided genome and phenotype can be achieved by used grammar)
    starting_inds = load_json_data()
    ge_instance.starting_inds = starting_inds

    pop = ge_instance.initialize_population(size=10, starting_inds=starting_inds)

    # valid starting individuals will remain in the properties of ge_instance
    assert ge_instance.starting_inds is not None
    assert len(pop) == 10


## Test for population initialization (with invalid starting ind)
def test_init_pop_invalid_starter_ind(ge_instance: GrammaticalEvolution):
    # initialize population with a invalid individual => initialize_population() should ignore this one and generate a random ind in its place
    starting_genome = [112, 14, 7, 10, 10, 5, 3, 2, 4, 1, 66, 55, 44]
    starting_invalid_phenotype = "AT(I_DONT_EXIST)-ACP--KS-AT(Aminomalonyl)-KR-DH-ACP--KS-AT(Hydroxymalonyl)-KR-DH-ER-ACP--KS-AT(Methoxymalonyl)-ACP--TE(Macrolactamization)"
    invalid_starting_ind = [starting_genome, starting_invalid_phenotype]
    ge_instance.starting_inds = [invalid_starting_ind]

    pop = ge_instance.initialize_population(size=10, starting_inds=[invalid_starting_ind])

    # the previously set starting ind should now have been removed, since it was not valid
    assert ge_instance.starting_inds is None
    assert len(pop) == 10


## Test for total Individual counter
def test_counter_total_ind(ge_instance: GrammaticalEvolution):
    # Inside the fitness evaluation the total Individual count (and also unique Individual count) is done, since 
    # that is one of the most expensive parts of this GE algorithm. Hence why we need to call evaluate_fitness() in this test case
    init_pop_size = 6
    # don't do something like this in the actual code, change the parameter in the config file instead!
    # We're doing it here since we need to ensure this number is odd, since the "parent generation" (aka onepoint_crossover function)
    # will always return an odd number and that could mess up our test case here
    evo_config.GENERATION_SIZE = 10

    pop_0 = ge_instance.initialize_population(size=init_pop_size, starting_inds=None)

    # trigger the fitness evaluation (which normally is done inside the evolutionary loop function)
    ge_instance.perform_evaluation(pop_0)

    # total Individual count should now be the same as the size of the initial population
    assert ge_instance.counter_total_ind == init_pop_size

    # now let's simulate one additional evolutionary step and recount the total Ind counter
    pop_0.sort(reverse=True) # we need to sort the initial population manually, since Gen 0 didn't run through the evolutionary loop
    pop_1, best_ever = ge_instance.evo_step(curr_pop=pop_0,
                                            best_ever=pop_0[0], 
                                            loop_count=1)

    assert ge_instance.counter_total_ind == init_pop_size + evo_config.GENERATION_SIZE


## Test for unique Individual counter (and counter for Individiual.evaluate() calls)
def test_unique_ind_counter(ge_instance: GrammaticalEvolution):
    init_pop_size = 6
    pop = ge_instance.initialize_population(size=init_pop_size, starting_inds=None)
    actual_unique_genomes = dict()
    valid_phenotypes = 0

    # Note: since individuals are created randomly, we need to check manually if all of them are actually unique
    for ind in pop:
        genome_tuple = tuple(ind.genome)  # Convert the list to a tuple, so that we can use it as index in the dictionary
        actual_unique_genomes[genome_tuple] = ind

    if len(actual_unique_genomes) == init_pop_size:
        # This should almost always be the case, since probability of identical Ind during initialization is rather small
        pop[1] = copy.copy(pop[0]) # let's add a duplicate by adding the first Ind also into the second position

    # Trigger fitness evaluation, so that the counter for unique and total Individuals are activated
    ge_instance.perform_evaluation(pop)

    # Loop over population again and check for which individuals we now have a phenotype
    for ind in pop:
        if ind.phenotype is not None:
            valid_phenotypes += 1
    
    # Since we have one duplicate, we have one less unique individual now
    assert len(ge_instance.unique_genomes) == init_pop_size - 1
    # Since we have one duplicate, we should have one less call (or possibly 2, if the duplicates have no valid phenotype) to Individual.evaluate()
    assert ge_instance.counter_eval_fit < valid_phenotypes

