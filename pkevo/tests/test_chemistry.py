import pytest
import os
import sys
import pprint

if sys.platform.startswith("win"):
    pytest.skip("skipping MOD-related tests on Windows", allow_module_level=True)

from chemistry.chemistry import Chemistry


@pytest.fixture
def chem_instance():
    chem = Chemistry(verbose=False)

    return chem


# Basic test to see if an instance of the Chemistry class can be initialized properly
def test_chem_init(chem_instance: Chemistry):
    """
    This test case evaluates the initialization process of the chemistry space, for which the MØD (https://github.com/jakobandersen/mod) 
    tool is used. It loads available molecules, reactions and substrategies. Note: No actual polyketides are created here in this step. 
    That is done in the next test case.
    """
    assert chem_instance is not None
    assert len(chem_instance.molecules['starter_units']) > 0
    assert len(chem_instance.molecules['extender_units']) > 0
    assert len(chem_instance.molecules['basic_units']) > 0
    assert len(chem_instance.reactions['fwd']) > 0
    assert len(chem_instance.reactions['rev']) == 0 # Note: Deconstructing molecules ATM not implemented, hence no rev. reactions
    assert len(chem_instance.substrategies['fwd']) > 0
    assert len(chem_instance.substrategies['rev']) == 0 # Note: Deconstructing molecules ATM not implemented hence no rev. subs

    # check if all molecule files were loaded correctly
    # from current path switch into the chemistry folder for molecules
    current_dir = os.path.dirname(os.path.abspath(__file__))

    molecules_folder = os.path.join(current_dir, '..', "chemistry", "chemistry_model", "molecules")

    molecule_files_count = 0
    # Iterate directory and subdirectories
    for root_dir, cur_dir, files in os.walk(molecules_folder):
        molecule_files_count += len(files)
    
    assert (len(chem_instance.molecules['starter_units']) + 
            len(chem_instance.molecules['extender_units']) + 
            len(chem_instance.molecules['basic_units'])) == molecule_files_count

    ### Note: Uncomment the following print statements if you're interested in seeing how the chemistry looks like
    ### Then run the pytest command with the additional '-s' argument to get the print statements displayed in the test results
    # print("Starter Units:")
    # print(chem_instance.molecules['starter_units'])

    # print("\nExtender Units:")
    # print(chem_instance.molecules['extender_units'])

    # print("\nBasic Units:")
    # print(chem_instance.molecules['basic_units'])

    # print("\nReactions:")
    # pprint.pprint(chem_instance.reactions['fwd'])
    # pprint.pprint(chem_instance.reactions['rev'])

    # print("\nSubstrategies")
    # pprint.pprint(chem_instance.substrategies)


def test_construct_molecule(chem_instance: Chemistry):
    """
    After the chemistry setup is complete, we can now try to create molecules by again using MØD. As input we use a Phenotype
    (more precisely the instruction notes for that Phenotype), for which we know that it was able to produce molecules in the past.
    In a first step, based on the available Substrategies, we now generate a 'complete' Strategy which is specific to that Phenotype.
    After that, we get MØD involved by asking it to create a derivation graph based on the molecules, reactions and that one
    specific strategy. In that graph, we should then be able to find our desired product(s). 
    """
    # Equivalent PKS Domain string: "AT(2-Methylbutyryl)-ACP--KS-AT(Aminomalonyl)-KR-DH-ER-ACP--KS-AT(Chlorethylmalonyl)-KR-DH-ACP--KS-AT(Ethylmalonyl)-ACP--KS-AT(Methoxymalonyl)-KR-DH-ER-ACP--TE(Hydrolysis)-Cyclization(3)"
    # To generate an assembly line directly out of a PKS domain string, we could also use `Grammar.deconstruct_string`. However, the Grammar class
    # is not involved in this test case, which should only be focusing on the Chemistry
    assembly_line = [
        'AT', '2-Methylbutyryl', 'ACP', 
        'KS', 'AT', 'Aminomalonyl', 'KR', 'DH', 'ER', 'ACP',
        'KS', 'AT', 'Chlorethylmalonyl', 'KR', 'DH', 'ACP',
        'KS', 'AT', 'Ethylmalonyl', 'ACP',
        'KS', 'AT', 'Methoxymalonyl', 'KR', 'DH', 'ER', 'ACP', 
        'TE', 'Hydrolysis', 'Cyclization', '3']
    
    expected_reactions_performed = [
        'Add 2-Methylbutyryl', 'Reduce chem. space', 'Add Aminomalonyl', 'Claisen condensation', 'Beta carbonyl reduction', 
        'Beta hydroxyl dehydration', 'Enoyl reduction', 'Reduce chem. space', 'Add Chlorethylmalonyl', 'Claisen condensation', 
        'Beta carbonyl reduction', 'Beta hydroxyl dehydration', 'Reduce chem. space', 'Add Ethylmalonyl', 'Claisen condensation', 
        'Reduce chem. space', 'Add Methoxymalonyl', 'Claisen condensation', 'Beta carbonyl reduction', 'Beta hydroxyl dehydration', 
        'Enoyl reduction', 'Reduce chem. space', 'Release reactions: Hydrolysis', 'Number of cyclizations: 3']

    # Create all possible Molecules (as SMILES strings) for the PKS assembly line instructions (given available molecules, reactions and strategies in Chemistry)
    result_smiles, performed_reactions, dg = chem_instance.get_smiles_for_ind(assembly_line, get_reactions=True)

    assert result_smiles is not None
    assert 'C(C(CC(C(C(CCCl)=CC(CC(C)CC)N)=O)CC)OC)(O)=O' in result_smiles

    assert performed_reactions == expected_reactions_performed


def test_triketide_lactone_generation(chem_instance: Chemistry):
    """
    Just like the previous test, but this time with a 'real-life' bio-engineered simple PKS, which should be capable of producing
    Triketide lactone. 
    Source: 2004 paper by Kira J. Weissman (https://pubmed.ncbi.nlm.nih.gov/15539364/)
    """
    # equivalent PKS Domain string: "AT(Propionyl)-ACP--KS-AT(Malonyl)-KR-ACP--KS-AT(Malonyl)-KR-ACP--TE(Macrolactonization)"
    assembly_line = [
        'AT', 'Propionyl', 'ACP',
        'KS', 'AT', 'Malonyl', 'KR', 'ACP',
        'KS', 'AT', 'Malonyl', 'KR', 'ACP', 
        'TE', 'Macrolactonization']
    
    expected_result_molecule = "C1(CC(CC(CC)O1)O)=O"

    expected_reactions_performed = [
        'Add Propionyl', 'Reduce chem. space', 'Add Malonyl', 'Claisen condensation', 'Beta carbonyl reduction', 
        'Reduce chem. space', 'Add Malonyl', 'Claisen condensation', 'Beta carbonyl reduction', 'Reduce chem. space', 
        'Release reactions: Macrolactonization']

    # Create all possible Molecules (as SMILES strings) for the PKS assembly line instructions (given available molecules, reactions and strategies in Chemistry)
    result_smiles, performed_reactions, dg = chem_instance.get_smiles_for_ind(assembly_line, get_reactions=True)

    assert result_smiles is not None
    assert expected_result_molecule in result_smiles
    assert performed_reactions == expected_reactions_performed


def test_geldanamycin_generation(chem_instance: Chemistry):
    # Eq. PKS Domainstring: "AT(AHBA)-ACP--KS-AT(Methylmalonyl)-KR-DH-ER-ACP--KS-AT(Methoxymalonyl)-KR-DH-ER-ACP--KS-AT(Methylmalonyl)-KR-ACP--KS-AT(Methylmalonyl)-KR-DH-ACP--KS-AT(Methoxymalonyl)-KR-ACP--KS-AT(Malonyl)-KR-DH-ER-ACP--KS-AT(Methylmalonyl)-KR-DH-ACP--TE(Macrolactamization)"
    assembly_line = [
        'AT', 'AHBA', 'ACP', 
        'KS', 'AT', 'Methylmalonyl', 'KR', 'DH', 'ER', 'ACP', 
        'KS', 'AT', 'Methoxymalonyl', 'KR', 'DH', 'ER', 'ACP', 
        'KS', 'AT', 'Methylmalonyl', 'KR', 'ACP', 
        'KS', 'AT', 'Methylmalonyl', 'KR', 'DH', 'ACP', 
        'KS', 'AT', 'Methoxymalonyl', 'KR', 'ACP', 
        'KS', 'AT', 'Malonyl', 'KR', 'DH', 'ER', 'ACP', 
        'KS', 'AT', 'Methylmalonyl', 'KR', 'DH', 'ACP', 
        'TE', 'Macrolactamization']
    
    expected_result_molecule = "C1CCC(C(C(C)=CC(C)C(C(CC(C)CC2:C:C(:C:C(:C:2)NC(C=1C)=O)O)OC)O)O)OC"
    
    expected_reactions_performed = [
        'Add AHBA', 'Reduce chem. space', 
        'Add Methylmalonyl', 'Claisen condensation', 'Beta carbonyl reduction', 'Beta hydroxyl dehydration', 'Enoyl reduction', 'Reduce chem. space', 
        'Add Methoxymalonyl', 'Claisen condensation', 'Beta carbonyl reduction', 'Beta hydroxyl dehydration', 'Enoyl reduction', 'Reduce chem. space', 
        'Add Methylmalonyl', 'Claisen condensation', 'Beta carbonyl reduction', 'Reduce chem. space', 
        'Add Methylmalonyl', 'Claisen condensation', 'Beta carbonyl reduction', 'Beta hydroxyl dehydration', 'Reduce chem. space', 
        'Add Methoxymalonyl', 'Claisen condensation', 'Beta carbonyl reduction', 'Reduce chem. space', 
        'Add Malonyl', 'Claisen condensation', 'Beta carbonyl reduction', 'Beta hydroxyl dehydration', 'Enoyl reduction', 'Reduce chem. space', 
        'Add Methylmalonyl', 'Claisen condensation', 'Beta carbonyl reduction', 'Beta hydroxyl dehydration', 'Reduce chem. space', 
        'Release reactions: Macrolactamization']
    
    result_smiles, performed_reactions, dg = chem_instance.get_smiles_for_ind(assembly_line, get_reactions=True)

    assert result_smiles is not None
    assert expected_result_molecule in result_smiles
    assert performed_reactions == expected_reactions_performed


def test_erythromycin(chem_instance: Chemistry):
    erythromycin_precursor_6deb = 'C1(C(C)C(C(C)C(C(C)CC(C(C(C)C(C(C)C(CC)O1)O)=O)C)O)O)=O'
    # Represents PKS Domainstring: AT(Propionyl)-ACP--KS-AT(Methylmalonyl)-KR-ACP--KS-AT(Methylmalonyl)-KR-ACP--KS-AT(Methylmalonyl)-ACP--KS-AT(Methylmalonyl)-KR-DH-ER-ACP--KS-AT(Methylmalonyl)-KR-ACP--KS-AT(Methylmalonyl)-KR-ACP--TE(Macrolactonization)
    assembly_line = [
        'AT', 'Propionyl', 'ACP', 
        'KS', 'AT', 'Methylmalonyl', 'KR', 'ACP', 
        'KS', 'AT', 'Methylmalonyl', 'KR', 'ACP', 
        'KS', 'AT', 'Methylmalonyl', 'ACP', 
        'KS', 'AT', 'Methylmalonyl', 'KR', 'DH', 'ER', 'ACP', 
        'KS', 'AT', 'Methylmalonyl', 'KR', 'ACP', 
        'KS', 'AT', 'Methylmalonyl', 'KR', 'ACP', 
        'TE', 'Macrolactonization']
    
    expected_reactions_performed = [
        'Add Propionyl', 'Reduce chem. space', 'Add Methylmalonyl', 'Claisen condensation', 'Beta carbonyl reduction', 'Reduce chem. space', 
        'Add Methylmalonyl', 'Claisen condensation', 'Beta carbonyl reduction', 'Reduce chem. space', 
        'Add Methylmalonyl', 'Claisen condensation', 'Reduce chem. space', 
        'Add Methylmalonyl', 'Claisen condensation', 'Beta carbonyl reduction', 'Beta hydroxyl dehydration', 'Enoyl reduction', 'Reduce chem. space', 
        'Add Methylmalonyl', 'Claisen condensation', 'Beta carbonyl reduction', 'Reduce chem. space', 
        'Add Methylmalonyl', 'Claisen condensation', 'Beta carbonyl reduction', 'Reduce chem. space', 
        'Release reactions: Macrolactonization']
    
    # Create all possible Molecules (as SMILES strings) for the PKS assembly line instructions (given available molecules, reactions and strategies in Chemistry)
    result_smiles, performed_reactions, dg = chem_instance.get_smiles_for_ind(assembly_line, get_reactions=True)

    # print(result_smiles)
    # print(performed_reactions)
    # print(len(result_smiles))

    assert result_smiles is not None
    assert erythromycin_precursor_6deb in result_smiles
    assert performed_reactions == expected_reactions_performed


def print_dg(dg):
    import subprocess
    from external_tools.mod_importer import mod
    # print derivation graph
    dg.print()

    # flush summary file handle
    mod.post.flushCommands()
    # generate summary/summery.pdf
    subprocess.run(["/home/talax/xtof/local/Mod/bin/mod_post"]) 
