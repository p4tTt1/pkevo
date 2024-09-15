from grammatical_evolution.ge_utils import crossovers
from grammatical_evolution.model._individual import _Individual
from config.pkevo_config import EvolutionConfig as evo_config


def test_variable_onepoint_crossover():
    """
    Test the one-point crossover operation (which is essentially responsible for the reproduction of individuals).
    Note: By default the resulting genomes may differ in length, since for each parent the crossover point is selected
    individually.
    """
    ind1 = _Individual(optimization_type="max", default_fitness=0.0)
    ind2 = _Individual(optimization_type="max", default_fitness=0.0)

    # Ensure that there will be a crossover by tampering with the config parameters
    evo_config.CROSSOVER_PROBABILITY = 1.0
    child1, child2 = crossovers.onepoint_crossover(ind1, ind2, default_fitness=0.0, variable_len=True)

    # Check that we got 2 Individuals as result
    assert isinstance(child1, _Individual)
    assert isinstance(child2, _Individual)
    # Check that overall genome length between parents and children stayed the same
    assert (len(ind1.genome) + len(ind2.genome)) == (len(child1.genome) + len(child2.genome))
    # Since crossover mutation is assured, check that child1 has first element of ind1 and last element of ind2 (and vice versa for child2)
    assert (ind1.genome[0] == child1.genome[0]) and (ind2.genome[-1] == child1.genome[-1])
    assert (ind2.genome[0] == child2.genome[0]) and (ind1.genome[-1] == child2.genome[-1])


def test_fixed_onepoint_crossover():
    """
    Test the one-point crossover operation with a fixed genome lenght. Both parents will have the crossover point at the
    same index, which ensures that the genome lenght of the offspring will stay the same.
    """
    ind1 = _Individual(optimization_type="max", length=20, default_fitness=0.0)
    ind2 = _Individual(optimization_type="max", length=20, default_fitness=0.0)

    evo_config.CROSSOVER_PROBABILITY = 1.0
    child1, child2 = crossovers.onepoint_crossover(ind1, ind2, default_fitness=0.0, variable_len=False)

    # check that we got 2 _Individuals as result
    assert isinstance(child1, _Individual)
    assert isinstance(child2, _Individual)
    # each genome should have the same lenght
    assert len(ind1.genome) == len(ind2.genome) == len(child1.genome) == len(child2.genome)


def test_onepoint_crossover_0_prob():
    """
    Test the one-point crossover operation with a 0% chance of crossovers happening => basically leading to cloning
    """
    ind1 = _Individual(optimization_type="max", default_fitness=0.0)
    ind2 = _Individual(optimization_type="max", default_fitness=0.0)

    # ensure that no crossovers will occur
    evo_config.CROSSOVER_PROBABILITY = 0.0
    child1, child2 = crossovers.onepoint_crossover(ind1, ind2, 0.0)

    # check that we got 2 _Individuals as result
    assert isinstance(child1, _Individual)
    assert isinstance(child2, _Individual)
    # check that the genomes of the children are identical to those of the parents
    assert ind1.genome == child1.genome
    assert ind2.genome == child2.genome


def test_variable_twopoint_crossover():
    """
    Test the two-point crossover operation with a variable genome lenght. As the name suggests, this operation selects two points
    on the genomes (for each parents different points) and the middle parts of the genomes are swapped out.
    Example: Parent 1: xxx|xx|xxxxx (len=10) -> Offspring a: xxxyyyyyxxxxx (len=13)
             Parent 2: yyyyy|yyyyy|yyyyyyy (len=17) -> Offspring b: yyyyyxxyyyyyyy (len=14)
    """
    ind1 = _Individual(optimization_type="max", genome=[1,1,1,1,1,1,1,1,1,1], default_fitness=0.0)
    ind2 = _Individual(optimization_type="max", genome=[2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2], default_fitness=0.0)

    evo_config.CROSSOVER_PROBABILITY = 1.0
    child1, child2 = crossovers.twopoint_crossover(ind1, ind2, default_fitness=0.0, variable_len=True)

    # Combined length of parental genomes should be equal to length of combined offspring genomes
    assert (len(ind1.genome) + len(ind2.genome)) == (len(child1.genome) + len(child2.genome))
    # Since crossover was guaranteed to happen, genomes of parents and children should differ
    assert ind1.genome != child1.genome


def test_fixed_twopoint_crossover():
    """
    Test the two-point crossover operation with a fixed genome lenght. As the name suggests, this operation selects two points
    on the genomes (on both parents the same points are selected) and the middle parts of the genomes are swapped out.
    Example: Parent 1: xxx|xx|xxxxx (len=10) -> Offspring a: xxxyyxxxxx (len=10)
             Parent 2: yyy|yy|yyyyy (len=10) -> Offspring b: yyyxxyyyyy (len=10)
    """
    ind1 = _Individual(optimization_type="max", genome=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1], default_fitness=0.0)
    ind2 = _Individual(optimization_type="max", genome=[2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2], default_fitness=0.0)

    evo_config.CROSSOVER_PROBABILITY = 1.0
    child1, child2 = crossovers.twopoint_crossover(ind1, ind2, default_fitness=0.0, variable_len=False)

    # Length of child genome should be the same as that of the parent
    assert len(ind1.genome) == len(child1.genome)
    assert len(ind2.genome) == len(child2.genome)
    # Since crossover was guaranteed to happen, genomes of parents and children should differ
    # Note: there is however still a small chance, that the genomes won't differ. This would happen when the 2 randomly chosen 
    # crossover points are the same points on the genome. Then this assertion would fail, however the function is working as 
    # intended nonetheless. Simply re-run the test case to check again.
    assert ind1.genome != child1.genome
