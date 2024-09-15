import copy

from config.pkevo_config import EvolutionConfig as evo_config
import config.custom_logger as custom_logger

logger = custom_logger.get_logger(__name__, "ge_logs.log")


def generational_replacement(new_pop, old_pop):
    """
    The self.elite_size best individuals (sorted array of the previous generation) 
    are appended to the new population if they are better than the worst individuals in new
    population.

    Comparing to `steady_state_replacement()`: Here we return the new population, potentially extended
    with best X individuals from old population (if those are better than worst new individuals)

    Attributes:
        new_pop (List[_Individual]): List containing the individuals of the newly created popuplation
        old_pop (List[_Individual]): List containing the individuals of the previous generation    
    Returns: 
        List: list of selected individuals for current generation in sorted order (best -> worst)
    """
  
    for ind in old_pop[:evo_config.ELITE_SIZE]:
        new_pop.append(copy.copy(ind))
    new_pop.sort(reverse=True)
    return new_pop[:evo_config.GENERATION_SIZE]


def steady_state_replacement(new_pop, old_pop):
    """
    Replaces the worst old individual with the best new individual, regardless of whether this new
    individual is actually better then the other.

    Comparing to `generational_replacement()`: Here we return the old population, which now contains
    the best new individuals (which replaced the worst old individiual, regardless of which of those two
    are actually better)
    
    Attributes:
        new_pop (List[_Individual]): List containing the individuals of the newly created popuplation
        old_pop (List[_Individual]): List containing the individuals of the previous generation   
    
    Returns:
        List: list of selected individuals for current generation in sorted order (best -> worst)
    """

    old_pop[-1] = max(new_pop)
    old_pop.sort(reverse=True)
    return old_pop