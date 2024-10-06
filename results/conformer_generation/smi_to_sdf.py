# Credits to Maximilian Faissner for initial script. Based on CDPKit Python Cookbook for 
# "Conformer ensembles": https://cdpkit.org/v1.0.0/cdpl_python_cookbook/confgen/ensembles.html 

import argparse
import os
import time

# CDPKit
import CDPL.Chem as Chem
import CDPL.ConfGen as ConfGen
import CDPL.Pharm as Pharm
from cdpkit_utils import generateConformationEnsembles

def main():
    start_time = time.perf_counter()
    parser = argparse.ArgumentParser(description="SMILES Data Processor")
    parser.add_argument("input_file", help="Input .smi file")
    args = parser.parse_args()

    input_file = args.input_file
    output_sdf = 'python_cookbook.sdf'
    max_time = 900  # in seconds
    min_rmsd = 0.5
    e_window = 15.0
    max_confs = 100
    writer = Chem.MolecularGraphWriter(output_sdf)

    try:
        with open(input_file, "r") as smi_file:
            for line in smi_file:
                line = line.strip()
                if not line:
                    continue

                smiles, name = line.split(" ", 1)
                compound_name = name.strip()
                mol = Chem.parseSMILES(smiles)
                Chem.setName(mol, compound_name)
                print(f'Processing {compound_name} with {mol.numAtoms} atoms and {mol.numBonds} bonds')

                conf_gen = ConfGen.ConformerGenerator()
                conf_gen.settings.timeout = max_time * 1000
                conf_gen.settings.minRMSD = min_rmsd
                conf_gen.settings.energyWindow = e_window
                conf_gen.settings.maxNumOutputConformers = max_confs

                mol_start_time = time.perf_counter()
                status, num_confs = generateConformationEnsembles(mol, conf_gen)
                mol_end_time = time.perf_counter()
                mol_elapsed_time = mol_end_time - mol_start_time
                print(f"Elapsed time (seconds) for molecule {compound_name}: {mol_elapsed_time:.3f}")

                if num_confs > 0:
                    if not writer.write(mol):
                        print(f'Error: output of conformer ensemble for molecule {compound_name} failed')
                    print(f'Wrote {num_confs} conformers into {output_sdf}')

    except FileNotFoundError:
        print(f"Error: File '{input_file}' not found.")
    except Exception as e:
        print(f"Error: An error occurred - {str(e)}")

    writer.close()
    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    print(f"Total elapsed time (seconds): {elapsed_time:.3f}")

if __name__ == "__main__":
    main()
