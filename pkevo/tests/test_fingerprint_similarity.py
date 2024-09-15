import pytest
from typing import List

from config.pkevo_config import SearchConfig as search_config

from grammatical_evolution.ge_utils.fitness_functions import FitnessFunctions
from grammatical_evolution.model import _individual
from grammatical_evolution.model import _molecule as _mol


@pytest.fixture
def individuals():
    """
    For convenience create some individuals and molecules and assign those to the individuals.
    """
    # Create some Individuals
    ind1 = _individual._Individual(optimization_type="max", default_fitness=-1.0)
    ind2 = _individual._Individual(optimization_type="max", default_fitness=-1.0)

    # Let's now create some molecules (during actual GE evolution simulation, this would be done by MOD)
    exact_molecule_clone = _mol.Molecule(smiles_string="C(C(C)CC)(C(C(O)=O)N)=O", fitness=0.0) # exact same as Target
    slightly_altered_molecule = _mol.Molecule(smiles_string="C(C(C)C)(C(C(O)=O)N)", fitness=0.0) # slightly altered compared to Target
    random_molecule = _mol.Molecule(smiles_string="Cc1ccccc1", fitness=0.0) # random SMILES
    molecule_for_ind2 = _mol.Molecule(smiles_string="c1c(N(=O)=O)cccc1", fitness=0.0) # another random SMILES

    # Assign molecules to Individuals (in code done by Chemistry.get_smiles_for_ind() function)
    ind1.molecules = [slightly_altered_molecule, random_molecule, exact_molecule_clone]
    ind2.molecules = [slightly_altered_molecule, random_molecule, exact_molecule_clone, molecule_for_ind2]

    return [ind1, ind2]


def test_molecule_evaluation(individuals: List[_individual._Individual]):
    """
    Test for checking the general functionality of the fitness function "objective_function_fingerprint_similarity".
    It should evaluate all molecules of an individual based on fingerprint similarity to a given target molecule and
    return the highest scored fitness value. As a side-effect, it should alter the Individual.molecules list by replacing
    the default fitness of each molecule with the calculated fitness 
    """
    # Define Target SMILES (in pkevo_config.yp called "TARGET_SMILES") and initialize empty dictionary for caching
    cache_dict = dict()
    search_config.TARGET_SMILES = "C(C(C)CC)(C(C(O)=O)N)=O"

    # For this test we only need one Individual
    ind = individuals[0]


    # Let's make sure that the molecules are stored in the Individual in exactly the same order as they were initialized in the list
    assert len(ind.molecules) == 3
    assert ind.molecules[0].smiles_string == "C(C(C)C)(C(C(O)=O)N)"
    assert ind.molecules[1].smiles_string == "Cc1ccccc1"
    # Also each of the molecules fitness score should still be 0.0
    assert ind.molecules[0].fitness == ind.molecules[1].fitness == ind.molecules[2].fitness == 0.0

    # Now let's perform the molecular fingerprint similarity fitness evaluation (which uses RDKit)
    fingerprint_fitness = FitnessFunctions.FingerprintSimilarity()
    ind.fitness = fingerprint_fitness.objective_function(individual=ind, remembered_similarities=cache_dict)

    # Check to see if order of the molecules remained the same
    assert ind.molecules[0].smiles_string == "C(C(C)C)(C(C(O)=O)N)"
    assert ind.molecules[1].smiles_string == "Cc1ccccc1"
    assert ind.molecules[2].smiles_string == "C(C(C)CC)(C(C(O)=O)N)=O"

    # The returned fitness scores should be between 0 and 1 (0% - 100% similarity)
    assert 0.0 <= ind.molecules[0].fitness <= 1.0
    assert 0.0 <= ind.molecules[1].fitness <= 1.0
    assert 0.0 <= ind.molecules[2].fitness <= 1.0

    # Now check if the highest scored Molecule.fitness is equal to the overall fitness score of the Individiual
    max_fitness_in_molecules = 0.0
    for mol in ind.molecules:
        if mol.fitness > max_fitness_in_molecules:
            max_fitness_in_molecules = mol.fitness

    assert max_fitness_in_molecules == ind.fitness
    # Note: You might be wondering why we waste memory space by storing the the overall best fitness of an individual twice 
    # (in Ind.fitness and in Ind.molecules[0].fitness) and that's a valid question. Reason is that, depending on selected fitness function,
    # there might not be any molecules involved in the simulation. On the other hand, Ind.fitness is not case-specific, but rather a universal
    # feature of GE and guaranteed to be present in each and every execution of this application for each and every individual.


def test_molecule_caching(individuals: List[_individual._Individual]):
    """
    To improve performance, already known and evaluated molecules are stored in a cache, which is a dictionary containing
    only unique molecules (unique is the SMILES representation of the molecule). This test checks if that dictionary is 
    correctly created and used.
    """
    # Define Target SMILES (in pkevo_config.yp called "TARGET_SMILES") and initialize empty dictionary for caching
    cache_dict = dict()
    search_config.TARGET_SMILES = "C(C(C)CC)(C(C(O)=O)N)=O"

    ind1 = individuals[0]
    ind2 = individuals[1]

    # Evaluate all molecules in first individual
    fingerprint_fitness = FitnessFunctions.FingerprintSimilarity()
    ind1.fitness = fingerprint_fitness.objective_function(individual=ind1, remembered_similarities=cache_dict)

    # Size of cached molecules should be equal to number of molecules in first Individual
    assert len(cache_dict) == len(ind1.molecules)

    # Evaluate all molecules of second individual
    ind2.fitness = fingerprint_fitness.objective_function(individual=ind2, remembered_similarities=cache_dict)

    # Size of cached molecules must be smaller than combined size of ind1.molecules and ind2.molecules, since some molecules are duplicates
    assert len(cache_dict) < (len(ind1.molecules) + len(ind2.molecules))

    # Finally, check to see if a fitness score returned by the cached dictionary (which for ind2 will be the case for any of the 3 duplicated 
    # molecules) actually corresponds to the fitness score of said molecule
    fit_calc_by_rdkit = ind1.molecules[0].fitness # this is the fitness score for the 'slightly altered molecule', which was calculated by RDKit
    fit_stored_in_cache = cache_dict[ind1.molecules[0].smiles_string] # this is stored fitness of said molecule in the cached dictionary
    fit_assigned_by_cache = ind2.molecules[0].fitness # for the 2. individual the same molecule must have the same fitness (this time assigned by the cache)
    assert fit_calc_by_rdkit == fit_stored_in_cache == fit_assigned_by_cache