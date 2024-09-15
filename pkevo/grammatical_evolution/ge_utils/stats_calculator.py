# Based on Robert Haas' "panakeias_garden" package (/panakeias_garden/pre_release_demos/metaheuristics.py)

import numpy as np
from collections import Counter

def _extract_basic_stats(values):
    """
    Calculate basic statistics from a list of values.

    Args:
        values (list): A list of numerical values.

    Returns:
        tuple: A tuple containing minimum, maximum, mean, standard deviation,
            10th percentile, 25th percentile, median, 75th percentile, and 90th percentile of the input values.
    """
    arr = np.array(values)
    return (
        np.min(arr),
        np.max(arr),
        np.mean(arr),
        np.std(arr),
        np.percentile(arr, 10),
        np.percentile(arr, 25),
        np.percentile(arr, 50),
        np.percentile(arr, 75),
        np.percentile(arr, 90)
    )


def calculate_statistics(individuals):
    """
    Calculate statistics for fitness and used codon values from a list of individuals.

    Args:
        individuals (list): A list of individual objects, each with fitness and used_codons attributes.

    Returns:
        tuple: 
            Two tuples containing statistics for fitness and used codon values.
                Each tuple contains minimum, maximum, mean, standard deviation,
                10th percentile, 25th percentile, median (50th percentile),
                75th percentile, and 90th percentile.
            Tuple containing number of unique, duplicated and invalid phenotypes
    """
    fitness_vals = []
    used_codon_vals = []

    # We are only interested in valid individuals, meaning those that expressed a phenotype
    valid_inds = [i for i in individuals if i.phenotype is not None]

    for i in valid_inds:
        fitness_vals.append(i.fitness)
        used_codon_vals.append(i.used_codons)
    
    fitness_stats = _extract_basic_stats(fitness_vals)
    used_codon_stats = _extract_basic_stats(used_codon_vals)
    
    # Calculate statistics for phenotypes
    valid_phenotypes = [i.phenotype for i in valid_inds] # valid phenotypes are all those that are not `None` 
    invalid_phenotypes = len(individuals) - len(valid_inds) # number of invalid phenotypes is equal to number of invalid inds
    unique_phenotypes = len(Counter(valid_phenotypes)) # number of truly unique phenotypes in provided list
    # Duplicated phenotypes are all those that are valid but not unique
    duplicated_phenotypes = len(valid_phenotypes) - unique_phenotypes
    phenotype_stats = (unique_phenotypes, duplicated_phenotypes, invalid_phenotypes)
        
    return fitness_stats, used_codon_stats, phenotype_stats
