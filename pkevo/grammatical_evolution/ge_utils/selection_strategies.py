import random

from config.pkevo_config import EvolutionConfig as evo_config


def tournament_selection(population, tournament_size=3):
    """
    Given an entire population, draw `tournament_size` competitors
    randomly and return the best.

    Arguments:
        population (List[_Individuals]): The current individuals of the population
        tournament_size (int, optional): Number defining how many individuals are drawn for each competition 

    Returns:
        List: List of the selected individuals
    """

    winners = []
    while len(winners) < evo_config.GENERATION_SIZE:
        competitors = random.sample(population, tournament_size)
        competitors.sort(reverse=True)
        winners.append(competitors[0])
    return winners


def truncation_selection(population, proportion=0.5):
    """
    Given an entire population, return the best `proportion` of
    them.

    Arguments:
        population (List[_Individuals]): The current individuals of the population
        proportion (float, optional): Float defining the share of how many individuals are selected

    Returns:
        List: List of the selected individuals    
    """

    population.sort(reverse=True)
    cutoff = int(len(population) * float(proportion))
    return population[0:cutoff]