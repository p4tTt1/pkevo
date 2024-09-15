from grammatical_evolution.model import _individual
from grammatical_evolution.model import _pks_architecture
import random


def test_individual_creation():
    # Create an _Individual instance with default parameters
    individual = _individual._Individual(optimization_type="max", default_fitness=0.0)

    # Test that attributes are initialized correctly
    assert isinstance(individual.genome, list)
    assert isinstance(individual.optimization_type, str)
    assert isinstance(individual.fitness, float)
    assert individual.phenotype is None
    assert len(individual.molecules) == 0
    assert isinstance(individual.used_codons, int)
    assert individual.assembly_line is None
    assert isinstance(individual.pet_name, str)

    # Create another individual with given properties
    individual2 = _individual._Individual(optimization_type="max", default_fitness=0.0, genome=[0,1,2,3], name_restriction="Something")

    assert individual2.genome == [0,1,2,3]
    assert "Something" in individual2.pet_name


def test_individual_fitness_comparison():
    # Test that for a optimization_type="max" setting, the individual with the higher fitness value is selected
    ind_max_1 = _individual._Individual(optimization_type="max", default_fitness=10.0)
    ind_max_2 = _individual._Individual(optimization_type="max", default_fitness=1.0)

    # Test that for a optimization_type="max" setting, the individual with the lower fitness value is selected
    ind_min_1 = _individual._Individual(optimization_type="min", default_fitness=10.0)
    ind_min_2 = _individual._Individual(optimization_type="min", default_fitness=1.0)

    assert ind_max_1 > ind_max_2
    assert max(ind_max_1, ind_max_2).fitness == 10
    assert ind_min_1 < ind_min_2
    assert max(ind_min_1, ind_min_2).fitness == 1


def test_population_sorting():
    # instead of performing actual fitness evaluation just use different default fitness values

    # first, let's test for mock fitness function with optimization type = min
    ind_min_1 = _individual._Individual(optimization_type="min", default_fitness=5.0)
    ind_min_2 = _individual._Individual(optimization_type="min", default_fitness=1.0)
    ind_min_3 = _individual._Individual(optimization_type="min", default_fitness=7.0)
    ind_min_4 = _individual._Individual(optimization_type="min", default_fitness=2.0)
    ind_min_5 = _individual._Individual(optimization_type="min", default_fitness=4.0)
    # merge individuals into a population
    min_pop = [ind_min_1, ind_min_2, ind_min_3, ind_min_4, ind_min_5]
    # sort() will call the overwritten __lt__() function from Individual
    # Note: usage of "reverse=True" is necessary due to the implemented logic in __lt__(), which in turn was chosen so that max() and min() will
    # behave in the expected way (those also rely on the overwritten __lt__() function) 
    min_pop.sort(reverse=True)

    # our list should now be sorted in ascending order with lowest Fitness at top (best to worst, where lower fitness score is better)
    assert min_pop[0].fitness == 1.0
    assert min_pop[1].fitness == 2.0
    assert min_pop[2].fitness == 4.0
    assert min_pop[3].fitness == 5.0
    assert min_pop[4].fitness == 7.0
    # just for sake of completion also ensure that the Individual objects are the same (ind_min_1 should be at 4th position in the population list now)
    assert ind_min_1 == min_pop[3]

    # now repeat same thing for the case where we use optimization type = max (since it affects how the __lt__() function works)
    ind_max_1 = _individual._Individual(optimization_type="max", default_fitness=5.0)
    ind_max_2 = _individual._Individual(optimization_type="max", default_fitness=1.0)
    ind_max_3 = _individual._Individual(optimization_type="max", default_fitness=7.0)
    max_pop = [ind_max_1, ind_max_2, ind_max_3]
    max_pop.sort(reverse=True)

    # our list should now be sorted in descending order with highest Fitness at top (best to worst, where higher fitness score is better)
    assert max_pop[0].fitness == 7.0
    assert max_pop[1].fitness == 5.0
    assert max_pop[2].fitness == 1.0
    assert ind_max_1 == max_pop[1]


def test_pks_architecture():
    # The PKS architecture is a more visual representation of the more commonly used PKS Domain string. 
    # It only serves the user to get a better understanding of the PKS structure necessary to build the desired polyketide.
    # The PKS architecture is not to be confused with the {assembly_line} property, which is a more verbose representation
    # (containing the rules used to achieve the phenotype), which is necessary for the Chemistry part of this application!

    # setup a bare minimum individual (its genotype is not relevant for this test case)
    ind = _individual._Individual(optimization_type="max", default_fitness=0.0, length=20)
    
    # first try to build an architecture directly (should not work since no phenotype exists yet)
    ind.build_pks_architecture()
    assert ind.architecture is None
    
    # give individual predifined phenotype (completly ignoring the actual genotype) 
    # consisting of 5 modules and 17 domains
    phenotype_string = "AT(AHBA)-ACP--KS-AT(Malonyl)-KR-DH-ACP--KS-AT(Malonyl)-KR-DH-ER-ACP--KS-AT(Methoxymalonyl)-ACP--TE(Claisen)"
    ind.phenotype = phenotype_string
    ind.build_pks_architecture()
    #print(ind.architecture.structure)

    assert isinstance(ind.architecture, _pks_architecture.PKSArchitecture)
    assert ind.architecture.module_count == 5
    assert ind.architecture.domain_count == 17
    assert ind.architecture.structure[0][1] == ['AT(AHBA)', 'ACP']
    assert ind.architecture.structure[1][1] == ['KS', 'AT(Malonyl)', 'KR', 'DH', 'ACP']
    assert ind.architecture.structure[2][1] == ['KS', 'AT(Malonyl)', 'KR', 'DH', 'ER', 'ACP'] 
    assert ind.architecture.structure[3][1] == ['KS', 'AT(Methoxymalonyl)', 'ACP'] 
    assert ind.architecture.structure[4][1] == ['TE(Claisen)']

    # finally let's try to "re-generate" the phenotype string by going through the pks_structure list sequentially
    # create a temporary copy of the original string without the dashes, those are not necessary
    tmp_pheno_str = phenotype_string.replace('-','')

    # sequentially, for each element in the list of lists, compare the beginning of the string with the current 
    # element and delete it from the string if it maches. In the end (if everything went well) we should have an empty string
    for module in ind.architecture.structure:
        for domain in module[1]: # Note: module[0] only contains a name for the module (e.g. "loading module", "1. module", ...)
            domain_name_len = len(domain)
            if domain == tmp_pheno_str[:domain_name_len]:
                tmp_pheno_str = tmp_pheno_str[domain_name_len:]

    assert tmp_pheno_str == ""

    

