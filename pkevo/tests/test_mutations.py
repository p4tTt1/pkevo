from grammatical_evolution.ge_utils import mutations
from grammatical_evolution.model._individual import _Individual
from config.pkevo_config import EvolutionConfig as evo_config


def test_singlepoint_mutation():
    # again ensure that we will have mutations by tampering with the config parameters
    evo_config.MUTATION_PROBABILITY = 1.0

    ind = _Individual(optimization_type="max", default_fitness=0.0)
    # create copy of original genome for later comparison
    original_genome = ind.genome.copy()

    # perfom single point mutations on individual
    ind = mutations.int_flip_mutation(ind)

    # since mutation probability is at 100%, the genomes should now differ from one another
    assert ind.genome != original_genome


def test_singlepoint_mutation_within_coding_reagion():
    # Again ensure that we will have mutations by tampering with the config parameters
    evo_config.MUTATION_PROBABILITY = 1.0

    ind = _Individual(optimization_type="max", length=80, default_fitness=0.0)

    # Let's assume that the first 20 genes of the genome are coding for a phenotype
    ind.used_codons = 20

    # Create copy of original genome for later comparison
    original_genome = ind.genome.copy()

    # Perfom single point mutations on individual
    ind = mutations.int_flip_mutation(ind)

    # Overall genome length was set to 80, first 20 genes should be different now,
    # whereas the remaining genome should look the same
    assert ind.genome[:20] != original_genome[:20]
    assert ind.genome[20:] == original_genome[20:]