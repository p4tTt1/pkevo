import random
import re

from config.pkevo_config import EvolutionConfig as evo_config
from grammatical_evolution.model._individual import _Individual


def onepoint_crossover(p_0: _Individual, p_1: _Individual, default_fitness, within_used=True, variable_len=True):
    """
    Perform one-point crossover between two parent individuals and return two children individuals. As the name suggests, the 
    crossover occurs at exactly one point of the individual's genome. How this point is selected depends on the arguments `within_used`
    and `variable_len`.

    Args:
        p_0 (_Individual): The first parent individual.
        p_1 (_Individual): The second parent individual.
        default_fitness: The default fitness value to be used for the offsprings.
        within_used (bool, optional): If True, selected crossover points will be within the used section of the genome. This ensures
                                      that the crossover will actually lead to changes, since it's affecting the coding region of the genome
                                      Defaults to True.
        variable_len (bool, optional): If True, for each parent individual the crossover point will be set separately. This can lead to 
                                       growing and shrinking genomes. If set to False, both parents will use the same crossover point, 
                                       which means that the genome lengths of the children will remain the same. Defaults to True.

    Returns:
        list: A list containing two offspring individuals resulting from one-point crossover.
    """
    # offspring creature should be of same creature-type as that of one of the parents (not a functional task, purely cosmetic purpose!)
    creatures = [re.search(r'-(\w+)$', name).group(1) for name in [p_0.pet_name, p_1.pet_name]]

    # Get the chromosomes
    genome_0, genome_1 = p_0.genome, p_1.genome

    # Uniformly generate crossover points. If within_used==True, selected points will be within the used section of the genomes.
    if within_used:
        max_p_0, max_p_1 = p_0.used_codons, p_1.used_codons
        # However, in case that genome has no used codons at all, default back to the overall length of the genome
        if max_p_0 == 0:
            max_p_0 = len(genome_0)
        if max_p_1 == 0:
            max_p_1 = len(genome_1)
    else:
        max_p_0, max_p_1 = len(genome_0), len(genome_1)

    # If genome length should be variable, we select a random point for each individual, which can lead to genomes growing or shrinking in size 
    if variable_len:
        # Example: genome 1: xxx|xxxxx (len=8)  ->  child a: xxxyyyy (len=7)
        #          genome 2: yyyy|yyyy (len=8)  ->  child b: yyyyxxxxx (len=9)
        pt_p_0, pt_p_1 = random.randint(1, max_p_0), random.randint(1, max_p_1)
    else:
        # same crossover point for both parents -> genome lenghts will remain the same
        pt_p_0 = pt_p_1 = random.randint(1, min(max_p_0, max_p_1)) # use min() function for precaution to avoid 'array idx out of bounds' errors 

    # Make new chromosomes by crossover only with the defined crossover probability
    if random.random() < evo_config.CROSSOVER_PROBABILITY:
        c_0 = genome_0[:pt_p_0] + genome_1[pt_p_1:]
        c_1 = genome_1[:pt_p_1] + genome_0[pt_p_0:]
    else:
        c_0, c_1 = genome_0[:], genome_1[:]

    # Put the new chromosomes into new individuals
    return [_Individual(p_0.optimization_type, default_fitness, genome=c_0, codon_size=evo_config.CODON_SIZE, name_restriction=random.choice(creatures)),
            _Individual(p_0.optimization_type, default_fitness, genome=c_1, codon_size=evo_config.CODON_SIZE, name_restriction=random.choice(creatures))]


def twopoint_crossover(p_0: _Individual, p_1: _Individual, default_fitness, within_used=True, variable_len=True):
    """
    Perform two-point crossover between two parent individuals and return two children individuals. As the name suggests, the 
    crossovers occur at two points of the individual's genome. How this point is selected depends on the arguments `within_used`.

    Args:
        p_0 (_Individual): The first parent individual.
        p_1 (_Individual): The second parent individual.
        default_fitness: The default fitness value to be used for the offsprings.
        within_used (bool, optional): If True, selected crossover points will be within the used section of the genome. This ensures
                                      that the crossover will actually lead to changes, since it's affecting the coding region of the genome
                                      Defaults to True.
        variable_len (bool, optional): If True, for each parent individual the crossover points will be set separately. This can lead to 
                                       growing and shrinking genomes in the offspring. If set to False, both parents will use the same 
                                       crossover point, which means that the genome lengths of the children will remain the same. 
                                       Defaults to True.

    Returns:
        list: A list containing two offspring individuals resulting from one-point crossover.
    """
    # offspring creature should be of same creature-type as that of one of the parents (not a functional task, purely cosmetic purpose!)
    creatures = [re.search(r'-(\w+)$', name).group(1) for name in [p_0.pet_name, p_1.pet_name]]

    # Get the parental genomes
    genome_0, genome_1 = p_0.genome, p_1.genome

    # Uniformly generate crossover points. If within_used==True, selected points will be within the used section of the genomes.
    if within_used:
        max_p_0, max_p_1 = p_0.used_codons, p_1.used_codons
        # However, in case that a genome has no used codons at all, default back to the overall length of the genome
        if max_p_0 == 0:
            max_p_0 = len(genome_0)
        if max_p_1 == 0:
            max_p_1 = len(genome_1)
    else:
        max_p_0, max_p_1 = len(genome_0), len(genome_1)

    # Each parent will have their own 2 crossing points
    if variable_len:
        first_cross, second_cross = random.randint(1, max_p_0), random.randint(1, max_p_0)
        third_cross, forth_cross = random.randint(1, max_p_1), random.randint(1, max_p_1)

        # Should by chance the two crossover points on the genome be the same, role the dice once more and than stick to the results
        # Note: Having the same crossover point twice basically means that a onepoint crossover is performed, which is not tragic in any way.
        # However, it is slightly inconvenient for the `test_crossover.py` test case, since there it will throw a false alert.
        if first_cross == second_cross:
            first_cross = random.randint(1, max_p_0)
        if third_cross == forth_cross:
            third_cross = random.randint(1, max_p_1)

        # Since we randomly drew two numbers, it might be that the first crossover point is actually higher then the other -> switch them out
        if first_cross > second_cross:
            tmp_var = first_cross
            first_cross = second_cross
            second_cross = tmp_var
        # Same check necessary for crossing points 3 and 4
        if third_cross > forth_cross:
            tmp_var = third_cross
            third_cross = forth_cross
            forth_cross = tmp_var
    else:
        # Same crossover points for both parents -> we need to pick lower of the 2 max values to avoid 'array idx out of bounds' errors
        upper_bound = min(max_p_0, max_p_1)
        first_cross, second_cross = random.randint(1, upper_bound), random.randint(1, upper_bound)
        
        # Again check the correct order of the crossing points
        if first_cross > second_cross:
            tmp_var = first_cross
            first_cross = second_cross
            second_cross = tmp_var
        # since we are using the fixed crossover variation, the second individual shall use the same points as the first individual   
        third_cross = first_cross
        forth_cross = second_cross

    # Make new chromosomes by crossover only with the defined crossover probability
    if random.random() < evo_config.CROSSOVER_PROBABILITY:
        c_0 = genome_0[:first_cross] + genome_1[third_cross:forth_cross] + genome_0[second_cross:]
        c_1 = genome_1[:third_cross] + genome_0[first_cross:second_cross] + genome_1[forth_cross:]
    else:
        c_0, c_1 = genome_0[:], genome_1[:]

    # Put the new chromosomes into new individuals
    return [_Individual(p_0.optimization_type, default_fitness, genome=c_0, codon_size=evo_config.CODON_SIZE, name_restriction=random.choice(creatures)),
            _Individual(p_0.optimization_type, default_fitness, genome=c_1, codon_size=evo_config.CODON_SIZE, name_restriction=random.choice(creatures))]
