from grammatical_evolution.model._molecule import Molecule

# a Molecule object in the GE context only contains a SMILES string representation (immutable property!) of said molecule and a fitness value, 
# which can get updated (usually the Molecule will be initialized prior to evaluation, so fitness will be set according
# to the default fitness of an Individual)
def test_init_and_update():
    # that's how a molecule will typically be initialized, with no fitness since evaluation has yet to occur
    mol = Molecule(smiles_string='C(=O)', fitness=0.0)

    # check if initialization was successfull
    assert isinstance(mol, Molecule)
    assert mol.fitness == 0.0
    assert mol.smiles_string == 'C(=O)'

    # let's assume the evaluation took place, so let's update the  Molecule object and check if it was successfull
    mol.fitness = 99.9
    assert mol.fitness == 99.9

    # smiles_string should be immutable (so that it can be used as a key in a dictionary). Reason: canonical SMILES are unique and
    # there is no need for us to work with the same molecule multiple times (one molecule will always get the same fitness score during a GE simulation)
    try:
        # let's pretend that by mistake someone tries to change the smiles_string. This should through an AttributeError
        # Note: of course one can always do: mol._smiles_string = "something"
        # to change the value, but there's nothing I can do about that. That's just how Python works...
        mol.smiles_string = 'I shall not be stored in this variable'
    except AttributeError:
        pass

    # trying to change the SMILES string via mol.smiles_string should not change the original value of that variable    
    assert mol.smiles_string == 'C(=O)'