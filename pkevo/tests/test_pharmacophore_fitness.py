import sys
import pytest

if sys.platform.startswith("win"):
    pytest.skip("skipping LigandScout-related tests on Windows", allow_module_level=True)

from grammatical_evolution.ge_utils.fitness_functions import FitnessFunctions
from grammatical_evolution.model import _individual
from grammatical_evolution.model import _molecule as _mol


@pytest.mark.long_running
def test_pharmacophore_fitness_function():
    """
    Test the overall functionality of the pharmacophore fitness function, which uses LigandScout (more precisely its 
    CLI tools `idbgen` to create a database of conformers from the passed on molecules and `iscreen` to perform the
    virtual screen of the target pharmacophore file against the created conformer-database).
    
    Note: This test case can take quite a while to finish (around 30-40 seconds).
    """
    # Create a couple of dummy individuals
    ind1 = _individual._Individual(optimization_type='max', default_fitness=-1.0)
    ind2 = _individual._Individual(optimization_type='max', default_fitness=-1.0)
    ind3 = _individual._Individual(optimization_type='max', default_fitness=-1.0)
    ind4 = _individual._Individual(optimization_type='max', default_fitness=-1.0)
    ind1.pet_name ="dog-1"
    ind2.pet_name ="cat-2"
    ind3.pet_name ="bird-3"
    ind4.pet_name ="cow-4"

    # Create and assign molecules to those individuals (consider multiple cases: 1 Ind -> n or zero Molecules, repeating Molecules)
    ind1_mols = [_mol.Molecule('C1(C(CCCl)C(C(C)C)O1)=O', 0.0),
                 _mol.Molecule('CC1CC(C(C(C=C(C(C(C=CC=C(C(=O)NC2=CC(=O)C(=C(C1)C2=O)OC)C)OC)OC(=O)N)C)C)O)OC', 0.0), # this should be Geldanamycin
                 _mol.Molecule('C(C(C(C(C(C(CCCl)C1C(C)C(C(C(C(C=C(C(O1)=O)N)OC)=O)O)O)=O)O)O)CCCl)(CC(O)=O)=O', 0.0)]
    ind3_mols = [ind1_mols[0],
                 _mol.Molecule('C(C(C(C(C(C(CCCl)C(C(C)C(C1C(C(C=C(C(O1)=O)N)OC)=O)O)O)=O)O)O)CCCl)(CC(O)=O)=O', 0.0)] 
    ind4_mols = [_mol.Molecule('C(CC(C(C(C1C(C(CCCl)C(C(C(C(C(C(C=C(C(O1)=O)N)OC)=O)O)O)C)O)=O)O)CCCl)=O)(O)=O', 0.0)]

    ind1.molecules = ind1_mols
    ind3.molecules = ind3_mols
    ind4.molecules = ind4_mols

    # Create a generation out of the Individuals (Note: Unlike other fitness functions, this one evaluates the generation as a whole instead of each
    # individual separately. Though strictly speaking, this has nothing to do with 'co-evolution' but is rather done for performance reasons)
    gen0 = [ind1, ind2, ind3, ind4]

    # As for all fitness functions, provide a cache dictionary for lookup of already known molecules (since target remains the same, a molecule
    # will always score the same, so no need to perform the docking simulation a second time).
    cache_dict = dict()

    # Perform the pharmacophore fitness evaluation on the generation
    pharma_fitness = FitnessFunctions.PharmacophoreScore()
    pharma_fitness.objective_function(individuals=gen0, remembered_pharma_scores=cache_dict)

    # Individuals with molecules should now have a fitness of >= 0,
    # Individuals with no molecules should remain at the defined default fitness of -1.0
    assert ind1.fitness >= 0
    assert ind2.fitness == -1.0
    assert ind3.fitness >= 0
    assert ind4.fitness >= 0

    # All molecules should have a fitness of 0.0 or higher, same molecules should
    #  have same fitness (the case for first molecules of Ind1 and Ind3)
    # Also the fitness of each individual should be equal to the highest fitness of its molecules
    for ind in gen0:
        ind_fitness = ind.fitness
        max_mol_fitness = 0.0 if len(ind.molecules) > 0 else -1.0 # init mol fitness with 0 if Ind has molecules, or -1.0 if Ind has no molecules
        for mol in ind.molecules:
            assert mol.fitness >= 0.0
            if mol.fitness > max_mol_fitness:
                max_mol_fitness = mol.fitness
        assert ind_fitness == max_mol_fitness


@pytest.mark.long_running
def Xtest_pharmacophore_fitness_function_per_individual():
    """
    This test's purpose is merely to demonstrate that it is not viable to perform a LigandScout simulation
    for each individual separately. This test is structured just like the one above, with the only difference being
    that we call the fitness function multiple times (each time with a "generation" consisting of exactly one individual).

    Note: If you want to check the performance for yourself, remove the 'X' from the function name (and rename above function, so 
    that it won't be executed) and run the test with pytest. 
    
    Note: This test case can take quite a while to finish (around 60 seconds).
    """
    # Create a couple of dummy individuals
    ind1 = _individual._Individual(optimization_type='max', default_fitness=-1.0)
    ind2 = _individual._Individual(optimization_type='max', default_fitness=-1.0)
    ind3 = _individual._Individual(optimization_type='max', default_fitness=-1.0)
    ind4 = _individual._Individual(optimization_type='max', default_fitness=-1.0)
    ind1.pet_name ="dog-1"
    ind2.pet_name ="cat-2"
    ind3.pet_name ="bird-3"
    ind4.pet_name ="cow-4"

    # Create and assign molecules to those individuals (consider multiple cases: 1 Ind -> n or zero Molecules, repeating Molecules)
    ind1_mols = [_mol.Molecule('C1(C(CCCl)C(C(C)C)O1)=O', 0.0),
                 _mol.Molecule('CC1CC(C(C(C=C(C(C(C=CC=C(C(=O)NC2=CC(=O)C(=C(C1)C2=O)OC)C)OC)OC(=O)N)C)C)O)OC', 0.0), # this should be Geldanamycin
                 _mol.Molecule('C(C(C(C(C(C(CCCl)C1C(C)C(C(C(C(C=C(C(O1)=O)N)OC)=O)O)O)=O)O)O)CCCl)(CC(O)=O)=O', 0.0)]
    ind3_mols = [ind1_mols[0],
                 _mol.Molecule('C(C(C(C(C(C(CCCl)C(C(C)C(C1C(C(C=C(C(O1)=O)N)OC)=O)O)O)=O)O)O)CCCl)(CC(O)=O)=O', 0.0)] 
    ind4_mols = [_mol.Molecule('C(CC(C(C(C1C(C(CCCl)C(C(C(C(C(C(C=C(C(O1)=O)N)OC)=O)O)O)C)O)=O)O)CCCl)=O)(O)=O', 0.0)]

    ind1.molecules = ind1_mols
    ind3.molecules = ind3_mols
    ind4.molecules = ind4_mols

    # We won't be using the generation for the actual fitness function call, but rather just for the evaluation afterwards
    gen0 = [ind1, ind2, ind3, ind4]

    # Define the binding target, against which our molecules will be tested
    pharmacophore_file = 'gdm_pharmacophore_relaxed.pmz' 
    # As for all fitness functions, provide a cache dictionary for lookup of already known molecules (since target remains the same, a molecule
    # will always score the same, so no need to perform the docking simulation a second time).
    cache_dict = dict()

    # Perform the pharmacophore fitness evaluation on the generation
    pharma_fitness = FitnessFunctions.PharmacophoreScore()
    pharma_fitness.objective_function(individuals=[ind1], pharmacophore_file=pharmacophore_file, remembered_pharma_scores=cache_dict)
    pharma_fitness.objective_function(individuals=[ind2], pharmacophore_file=pharmacophore_file, remembered_pharma_scores=cache_dict)
    pharma_fitness.objective_function(individuals=[ind3], pharmacophore_file=pharmacophore_file, remembered_pharma_scores=cache_dict)
    pharma_fitness.objective_function(individuals=[ind4], pharmacophore_file=pharmacophore_file, remembered_pharma_scores=cache_dict)

    # Individuals with molecules should now have a fitness of >= 0,
    # Individuals with no molecules should remain at the defined default fitness of -1.0
    assert ind1.fitness >= 0
    assert ind2.fitness == -1.0
    assert ind3.fitness >= 0
    assert ind4.fitness >= 0

    # All molecules should have a fitness of 0.0 or higher, same molecules should
    #  have same fitness (the case for first molecules of Ind1 and Ind3)
    # Also the fitness of each individual should be equal to the highest fitness of its molecules
    for ind in gen0:
        ind_fitness = ind.fitness
        max_mol_fitness = 0.0 if len(ind.molecules) > 0 else -1.0 # init mol fitness with 0 if Ind has molecules, or -1.0 if Ind has no molecules
        for mol in ind.molecules:
            assert mol.fitness >= 0.0
            if mol.fitness > max_mol_fitness:
                max_mol_fitness = mol.fitness
        assert ind_fitness == max_mol_fitness

