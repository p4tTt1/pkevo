# Taken from ModByExamples Tuturial created by Prof. Christoph Flamm (file: vizSMI.py)
# Execute like this: python -m scripts_for_analysis.visualize_mod_objects (from root_dir/pkevo/)
import subprocess
from rdkit import DataStructs
from rdkit import Chem
from rdkit.Chem import AllChem 

from external_tools.mod_importer import mod
from chemistry.chemistry import Chemistry

# chem = Chemistry(verbose=False)
def is_in_subset_or_in_basics(graph, state, is_first_call):
    """
        Predicate to check if a graph is in a subset or in basic molecules.
        Args:
            graph: The graph to check.
            state: The current state of the synthesis.
            is_first_call (bool): Indicates if it's the first call.
        Returns:
        bool: True if the graph is in the subset or in basic molecules, False otherwise.
    """
    return graph in state.subset or graph in basic_molecules

def remove_from_universe(molecule):
    return mod.filterUniverse(lambda g, gs, first_call: g != molecule)


propionyl = mod.smiles("[Ag]SC(=O)CC", name="Propionyl")
malonyl = mod.smiles("[Au]SC(=O)CC(=O)O", name="Malonyl")
water = mod.smiles("[H][O][H]", name="Water")
proton = mod.smiles("[H+]", name="Proton")
hydride = mod.smiles("[H-]", name="Hydride")
basic_molecules = [water, proton, hydride]

cl_cond = mod.ruleGMLString("""
rule [
  ruleID "Intermolecular decarboxylative Claisen thioester condensation"
  left [
    edge [ source 1 target 2 label "-" ]
    edge [ source 6 target 7 label "-" ]
    edge [ source 9 target 10 label "-" ]
    edge [ source 10 target 12 label "-" ]
    edge [ source 12 target 13 label "-" ]
  ]
  context [
    node [ id 0 label "Ag" ]
    node [ id 1 label "S" ]
    node [ id 2 label "C" ]
    node [ id 3 label "O" ]
    node [ id 4 label "C" ]
    edge [ source 0 target 1 label "-" ]
    edge [ source 2 target 3 label "=" ]
    edge [ source 2 target 4 label "-" ]
    node [ id 5 label "Au" ]
    node [ id 6 label "S" ]
    node [ id 7 label "C" ]
    node [ id 8 label "O" ]
    node [ id 9 label "C" ]
    node [ id 10 label "C" ]
    node [ id 11 label "O" ]
    node [ id 12 label "O" ]
    node [ id 13 label "H" ]
    edge [ source 5 target 6 label "-" ]
    edge [ source 7 target 8 label "=" ]
    edge [ source 7 target 9 label "-" ]
    edge [ source 10 target 11 label "=" ]
  ]
  right [
    edge [ source 9 target 2 label "-" ]
    edge [ source 1 target 7 label "-" ]
    edge [ source 6 target 13 label "-" ]
    edge [ source 10 target 12 label "=" ]
  ]
]

""")

carb_red = mod.ruleGMLString("""
rule [
  ruleID "Reduction of beta carbonyl group"
  left [
    node [ id 8 label "H-" ]
    node [ id 9 label "H+" ]
    edge [ source 5 target 6 label "=" ]
  ]
  context [
    node [ id 0 label "Ag" ]
    node [ id 1 label "S" ]
    node [ id 2 label "C" ]
    node [ id 3 label "O" ]
    node [ id 4 label "C" ]
    node [ id 5 label "C" ]
    node [ id 6 label "O" ]
    node [ id 7 label "C" ]
    edge [ source 0 target 1 label "-" ]
    edge [ source 1 target 2 label "-" ]
    edge [ source 2 target 3 label "=" ]
    edge [ source 2 target 4 label "-" ]
    edge [ source 4 target 5 label "-" ]
    edge [ source 5 target 7 label "-" ]
  ]
  right [
    node [ id 8 label "H" ]
    node [ id 9 label "H" ]
    edge [ source 5 target 6 label "-" ]
    edge [ source 5 target 8 label "-" ]
    edge [ source 6 target 9 label "-" ]
  ]
]

""")

hydrox_rel = mod.ruleGMLString("""
rule [
  ruleID "Release from the synthase by nucleophilic attack of a hydroxyl group"
  left [
    edge [ source 2 target 3 label "-" ]
    edge [ source 6 target 7 label "-" ]
  ]
  context [
    node [ id 1 label "Ag" ]
    node [ id 2 label "S" ]
    node [ id 3 label "C" ]
    node [ id 4 label "O" ]
    edge [ source 1 target 2 label "-" ]
    edge [ source 3 target 4 label "=" ]
    node [ id 6 label "O" ]
    node [ id 7 label "H" ]
  ]
  right [
    edge [ source 2 target 7 label "-" ]
    edge [ source 6 target 3 label "-" ]
  ]
]
""")


# # equivalent PKS Domain string: "AT(Propionyl)-ACP--KS-AT(Malonyl)-KR-ACP--KS-AT(Malonyl)-KR-ACP--TE(Macrolactonization)"
# assembly_line = [
#         'AT', 'Propionyl', 'ACP',
#         'KS', 'AT', 'Malonyl', 'KR', 'ACP',
#         'KS', 'AT', 'Malonyl', 'KR', 'ACP', 
#         'TE', 'Macrolactonization']
    
# expected_result_molecule = "C1(CC(CC(CC)O1)O)=O"

# expected_reactions_performed = [
#         'Add Propionyl', 'Reduce chem. space', 'Add Malonyl', 'Claisen condensation', 'Beta carbonyl reduction', 
#         'Reduce chem. space', 'Add Malonyl', 'Claisen condensation', 'Beta carbonyl reduction', 'Reduce chem. space', 
#         'Release reactions: Macrolactonization']

# # Create all possible Molecules (as SMILES strings) for the PKS assembly line instructions (given available molecules, reactions and strategies in Chemistry)
# result_smiles, performed_reactions, dg = chem.get_smiles_for_ind(assembly_line, get_reactions=True)
# print(result_smiles)

## switch to term rewite
ls = mod.LabelSettings(mod.LabelType.Term, mod.LabelRelation.Unification)
dg = mod.DG(graphDatabase=mod.inputGraphs, labelSettings=ls)

strat = (
    mod.addUniverse([water, proton, hydride]) 
    # AT(Propionyl)-ACP
    >> mod.addSubset(propionyl)
    # KS-AT(Malonyl)-KR-ACP 
    >> mod.filterUniverse(is_in_subset_or_in_basics) >> mod.addSubset(malonyl) >> cl_cond >> carb_red 
    # KS-AT(Malonyl)-KR-ACP 
    >> mod.filterUniverse(is_in_subset_or_in_basics) >> mod.addSubset(malonyl) >> cl_cond >> carb_red
    # TE(Macrolactonization) 
    >> mod.filterUniverse(is_in_subset_or_in_basics) >> mod.leftPredicate[lambda derivation: len(derivation.left) == 1](hydrox_rel))

dg.build().execute(strat)

# print derivation graph
dg.print()


# print input rule(s)
for r in mod.inputRules:
    p = mod.GraphPrinter()
    p.setReactionDefault()
    p.withIndex = True
    r.print(p)



# for m in mod.inputGraphs:
#     m.print()

# flush summary file handle
mod.post.flushCommands()
# generate summary/summery.pdf
subprocess.run(["/home/talax/xtof/local/Mod/bin/mod_post"])