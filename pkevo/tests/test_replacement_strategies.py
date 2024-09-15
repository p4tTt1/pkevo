from grammatical_evolution.ge_utils import replacement_strategies
from grammatical_evolution.model._individual import _Individual
from config.pkevo_config import EvolutionConfig as evo_config


def test_generational_replacement():
    """
    The generational replacement keeps the `ELITE_SIZE` best individuals from the previous generation,
    if they are fitter than the worst individuals from the new population.
    """
    # Create an initial population of individuals (aka Gen 0)
    init_pop = [
        _Individual('max', 0.6, None),
        _Individual('max', 0.8, None),
        _Individual('max', 0.5, None), # this one is the 4th best overall, however only 3rd best in this gen and should therefore be eliminated
        _Individual('max', 0.3, None),
        _Individual('max', 0.4, None),
    ]
    init_pop.sort(reverse=True) # necessary, since the sorting normally happens after the replacement step!

    # Create a new sample population (aka Gen 1)
    new_pop = [
        _Individual('max', 0.5, None),
        _Individual('max', 0.7, None),
        _Individual('max', 0.3, None),
        _Individual('max', 0.1, None),
        _Individual('max', 0.2, None),
    ]

    # This means that we expect the best 2 inds from old pop to make it into the new pop (at least if their fitness)
    evo_config.ELITE_SIZE = 2

    # Apply generational replacement
    evo_config.GENERATION_SIZE = 5
    evo_config.POPULATION_SIZE = 5
    new_pop_result = replacement_strategies.generational_replacement(new_pop, init_pop)

    # Test assertions
    assert len(new_pop_result) == evo_config.GENERATION_SIZE
    assert new_pop_result[0].fitness == 0.8 # old ind
    assert new_pop_result[1].fitness == 0.7 # new ind
    assert new_pop_result[2].fitness == 0.6 # old ind
    assert new_pop_result[3].fitness == 0.5 # new ind
    assert new_pop_result[4].fitness == 0.3 # new ind


def test_steady_state_replacement():
    """
    The steady-state replacement strategy should always replace the worst individual
    from the old generation with the best individual from the new population, no matter
    how their two fitness scores compare.
    """
    # Create an initial population of individuals (aka Gen 0)
    init_pop = [
        _Individual('max', 0.6, None),
        _Individual('max', 0.8, None),
        _Individual('max', 0.5, None),
        _Individual('max', 0.3, None),
        _Individual('max', 0.4, None),
    ]
    init_pop.sort(reverse=True) # necessary, since the sorting normally happens after the replacement step!

    # Create a new sample population (aka Gen 1)
    new_pop = [
        _Individual('max', 0.5, None), # this one will not make the final cut despite being better than some from the other generation
        _Individual('max', 0.7, None), # since we only consider THE best individual in this population
        _Individual('max', 0.3, None),
        _Individual('max', 0.1, None),
        _Individual('max', 0.2, None),
    ]

    # Apply steady-state replacement
    evo_config.GENERATION_SIZE = 5
    evo_config.POPULATION_SIZE = 5
    new_pop_result = replacement_strategies.steady_state_replacement(new_pop, init_pop)

    # Test assertions
    assert len(new_pop_result) == evo_config.GENERATION_SIZE
    assert new_pop_result[0].fitness == 0.8 # old ind
    assert new_pop_result[1].fitness == 0.7 # new ind
    assert new_pop_result[2].fitness == 0.6 # old ind
    assert new_pop_result[3].fitness == 0.5 # old ind
    assert new_pop_result[4].fitness == 0.4 # old ind
