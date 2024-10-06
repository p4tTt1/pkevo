# Credits to Maximilian Faissner for initial script.

import CDPL.Chem as Chem
import CDPL.ConfGen as ConfGen
import CDPL.Pharm as Pharm
import CDPL.Biomol as Biomol

import sys


def generate3dConformation(mol: Chem.Molecule, struct_gen: ConfGen.StructureGenerator) -> int:
    """
    Generates a low-energy 3D structure of the argument molecule using the provided initialized ConfGen.StructureGenerator instance.

    Parameters:
    - mol (Chem.Molecule): Molecule to generate a 3D structure for.
    - struct_gen (ConfGen.StructureGenerator): Instance of the ConfGen.StructureGenerator class.

    Returns:
    - int: Status code indicating the success of the 3D structure generation.
    """
    # prepare the molecule for 3D structure generation
    ConfGen.prepareForConformerGeneration(mol)

    # generate the 3D structure
    status = struct_gen.generate(mol)

    # if sucessful, store the generated conformer ensemble as
    # per atom 3D coordinates arrays (= the way conformers are represented in CDPKit)
    if status == ConfGen.ReturnCode.SUCCESS:
        struct_gen.setCoordinates(mol)

    # return status code
    return status


def generateConformationEnsembles(mol: Chem.BasicMolecule, conf_gen: ConfGen.ConformerGenerator) -> (int, int):
    """
    Generates a conformation ensemble for the argument molecule using the provided initialized ConfGen.ConformerGenerator instance.

    Parameters:
    - mol (Chem.BasicMolecule): Molecule to generate a conformation ensemble for.
    - conf_gen (ConfGen.ConformerGenerator): Instance of the ConfGen.ConformerGenerator class.

    Returns:
    - int: Status code indicating the success of the conformation ensemble generation.
    - int: Number of generated conformers.
    """
    # prepare the molecule for conformer generation
    ConfGen.prepareForConformerGeneration(mol)

    # generate the conformer ensemble
    status = conf_gen.generate(mol)
    num_confs = conf_gen.getNumConformers()

    # if sucessful, store the generated conformer ensemble as
    # per atom 3D coordinates arrays (= the way conformers are represented in CDPKit)
    if status == ConfGen.ReturnCode.SUCCESS or status == ConfGen.ReturnCode.TOO_MUCH_SYMMETRY:
        conf_gen.setConformers(mol)
    else:
        num_confs = 0

    # return status code and the number of generated conformers
    return (status, num_confs)


def generatePharmacophore(mol: Chem.Molecule) -> Pharm.Pharmacophore:
    """
    Generates the pharmacophore of the molecule.

    Parameters:
    - mol (Chem.Molecule): Molecule to generate a pharmacophore for.

    Returns:
    - Pharm.Pharmacophore: Pharmacophore of the argument molecule.
    """
    Pharm.prepareForPharmacophoreGeneration(
        mol)    # first call utility function preparing the molecule for pharmacophore generation

    # create an instance of the pharmacophore generator default implementation
    ph4_gen = Pharm.DefaultPharmacophoreGenerator()
    # create an instance of the default implementation of the Pharm.Pharmacophore interface
    ph4 = Pharm.BasicPharmacophore()
    # use the name of the input molecule as pharmacophore name
    ph4_name = Chem.getName(mol)

    ph4_gen.generate(mol, ph4)          # generate the pharmacophore
    Pharm.setName(ph4, ph4_name)        # set the pharmacophore name

    return ph4


# outputs all (available) properties of the features stored in the given feature container
def outputFeatureProps(ph4: Pharm.FeatureContainer) -> None:
    ftr_type_str = {Pharm.FeatureType.UNKNOWN: 'UNKNOWN',
                    Pharm.FeatureType.HYDROPHOBIC: 'HYDROPHOBIC',
                    Pharm.FeatureType.AROMATIC: 'AROMATIC',
                    Pharm.FeatureType.NEGATIVE_IONIZABLE: 'NEGATIVE_IONIZABLE',
                    Pharm.FeatureType.POSITIVE_IONIZABLE: 'POSITIVE_IONIZABLE',
                    Pharm.FeatureType.H_BOND_DONOR: 'H_BOND_DONOR',
                    Pharm.FeatureType.H_BOND_ACCEPTOR: 'H_BOND_ACCEPTOR',
                    Pharm.FeatureType.HALOGEN_BOND_DONOR: 'HALOGEN_BOND_DONOR',
                    Pharm.FeatureType.HALOGEN_BOND_ACCEPTOR: 'HALOGEN_BOND_ACCEPTOR',
                    Pharm.FeatureType.EXCLUSION_VOLUME: 'EXCLUSION_VOLUME'}

    geom_str = {Pharm.FeatureGeometry.UNDEF: 'UNDEF',
                Pharm.FeatureGeometry.SPHERE: 'SPHERE',
                Pharm.FeatureGeometry.VECTOR: 'VECTOR',
                Pharm.FeatureGeometry.PLANE: 'PLANE'}

    print('Composition of pharmacophore \'%s\':' % Pharm.getName(ph4))

    for i in range(0, len(ph4)):
        ftr = ph4[i]

        print(' - Feature #%s:' % str(i))
        print('  - Type: %s' % ftr_type_str[Pharm.getType(ftr)])
        print('  - Geometry: %s' % geom_str[Pharm.getGeometry(ftr)])
        print('  - Tolerance: %s' % Pharm.getTolerance(ftr))
        print('  - Weight: %s' % Pharm.getWeight(ftr))
        print('  - Optional: %s' % Pharm.getOptionalFlag(ftr))
        print('  - Disabled: %s' % Pharm.getDisabledFlag(ftr))
        print('  - Length: %s' % Pharm.getLength(ftr))
        print('  - Hydrophobicity: %s' % Pharm.getHydrophobicity(ftr))

        # Pharm.Feature derives from Chem.Entity3D - therefore a function from the Chem package is used here!
        if Chem.has3DCoordinates(ftr):
            print('  - Position: %s' % Chem.get3DCoordinates(ftr))

        if Pharm.hasOrientation(ftr):
            print('  - Orientation: %s' % Pharm.getOrientation(ftr))


def processReceptorStructure(path: str, strip_res_list: bool) -> Chem.Molecule:
    """
    Reads and preprocesses the specified receptor structure.

    Parameters:
    - path (str): Path to the receptor structure file.
    - strip_res_list (bool): Whitespace separated list of PDB three-letter codes specifying residues to remove from the receptor structure (e.g. an existing ligand).

    Returns:
    - Chem.Molecule: Receptor structure.

    """

    # create reader for receptor structure (format specified by file extension)
    reader = Chem.MoleculeReader("path_to_receptor_structure_file.pdb")

    sup_fmts = [ Chem.DataFormat.MOL2,
                Biomol.DataFormat.PDB,
                Biomol.DataFormat.MMTF ]

    # check if the format is supported by this script
    if reader.getDataFormat() not in sup_fmts:
        sys.exit('Error: receptor input file format \'%s\' not supported' % name_and_ext[1])

    rec_mol = Chem.BasicMolecule()    # create an instance of the default implementation of the
                                      # Chem.Molecule interface that will store the receptor struct.
    try:
        if not reader.read(rec_mol):  # read receptor structure
            sys.exit('Error: reading receptor structure failed')

    except Exception as e:
        sys.exit('Error: reading receptor structure failed:\n' + str(e))

    # preprocess the receptor structure (removal of residues and
    # calculation of properties required by the pharm. generation procedure)
    try:
        # if structure comes from an MOL2 file, convert MOL2 residue data into PDB-style data
        if reader.getDataFormat() == Chem.DataFormat.MOL2:
            Biomol.convertMOL2ToPDBResidueInfo(rec_mol, True)

        rem_atoms = False

        # delete atoms belonging to residues that should be stripped
        if strip_res_list:
            atoms_to_rem = Chem.Fragment() # will store the atoms to delete
            res_to_strip = { tlc.upper() for tlc in strip_res_list }

            for atom in rec_mol.atoms:     # identify and note atoms belonging to the stripped residues
                if Biomol.getResidueCode(atom).upper() in res_to_strip:
                    atoms_to_rem.addAtom(atom)

            if atoms_to_rem.numAtoms > 0:
                rec_mol -= atoms_to_rem    # delete atoms from the receptor structure
                rem_atoms = True

        # prepares the receptor structure for pharmacophore generation
        Chem.perceiveSSSR(rec_mol, rem_atoms)
        Chem.setRingFlags(rec_mol, rem_atoms)
        Chem.calcImplicitHydrogenCounts(rec_mol, rem_atoms)
        Chem.perceiveHybridizationStates(rec_mol, rem_atoms)
        Chem.setAromaticityFlags(rec_mol, rem_atoms)

        if Chem.makeHydrogenComplete(rec_mol):                    # make implicit hydrogens (if any) explicit
            Chem.calcHydrogen3DCoordinates(rec_mol)               # calculate 3D coordinates for the added expl. hydrogens
            Biomol.setHydrogenResidueSequenceInfo(rec_mol, False) # set residue information for the added expl. hydrogens

        MolProp.calcAtomHydrophobicities(rec_mol, False)          # calculate atom hydrophobicity values (needed for hydrophobic
                                                                # pharm. feature generation)
    except Exception as e:
        sys.exit('Error: processing of receptor structure failed: ' + str(e))

    return rec_mol
