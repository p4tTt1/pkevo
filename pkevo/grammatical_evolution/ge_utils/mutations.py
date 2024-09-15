import random

from config.pkevo_config import EvolutionConfig as evo_config
from grammatical_evolution.model._individual import _Individual


def int_flip_mutation(individual: _Individual, within_used=True):
    """
    Perform integer flip mutation on an individual's genome.

    Args:
        individual (_Individual): The individual to mutate.
        within_used (bool, optional): If True, mutations will only occur in the coding region of the genome, unless there are no coding genes, 
                                      in which case the whole genome is used. Defaults to True.

    Returns:
        _Individual: The mutated individual.
    """
    # Define the allowed region where mutations might later take place
    mutation_region = len(individual.genome) # mutation might occur on the whole genome of the individual

    # Mutating regions of the genome that do not encode for a phenotype (e.g. in the tail region of the genome) don't have any
    # immediate impact on the individual, so it can be a good idea to limit mutations to only the coding region of the genome
    if within_used and individual.used_codons > 0:
        mutation_region = individual.used_codons

    # Iterate over selected region of the genome codon by codon
    for i in range(mutation_region):
        # Mutate a codon based on given probability
        if random.random() < evo_config.MUTATION_PROBABILITY:
            individual.genome[i] = random.randint(0, evo_config.CODON_SIZE)

    return individual


def deletion_mutation(individual: _Individual, within_used=True):
    """
    Perform deletion mutation on an individual's genome.

    Args:
        individual (_Individual): The individual to mutate.
        within_used (bool, optional): If True, deletion will only occur in the coding region of the genome, unless there are no coding genes, 
                                      in which case the whole genome is used. Defaults to True.

    Returns:
        _Individual: The mutated individual.
    """
    # Define the allowed region where mutations might later take place
    mutation_region = len(individual.genome) # mutation might occur on the whole genome of the individual
    if within_used and individual.used_codons > 0:
        mutation_region = individual.used_codons

    # Determine the position of the codon to delete
    position_to_delete = random.randint(0, mutation_region - 1)
    
    # Delete the codon at the specified position
    del individual.genome[position_to_delete]
    
    return individual


def insertion_mutation(individual: _Individual, within_used=True):
    """
    Perform insertion mutation on an individual's genome.

    Args:
        individual (_Individual): The individual to mutate.
        within_used (bool, optional): If True, deletion will only occur in the coding region of the genome, unless there are no coding genes, 
                                      in which case the whole genome is used. Defaults to True.

    Returns:
        _Individual: The mutated individual.
    """
    # Define the allowed region where mutations might later take place
    mutation_region = len(individual.genome) # mutation might occur on the whole genome of the individual
    if within_used and individual.used_codons > 0:
        mutation_region = individual.used_codons

    # Determine the position for inserting a new codon
    position_to_insert = random.randint(0, mutation_region)
    
    # Generate a new random codon to insert
    new_codon = random.randint(0, evo_config.CODON_SIZE)

    # Insert the new codon at the specified position
    individual.genome.insert(position_to_insert, new_codon)

    return individual