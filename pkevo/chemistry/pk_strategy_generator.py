# Note: Based on Robert Haas' panakeias_garden package (file: _enzyme_chemistry.py)
from rdkit import Chem
from rdkit.Chem import Descriptors
from indigo import Indigo

from external_tools.mod_importer import mod as _mod
import config.custom_logger as custom_logger

from .chemistry_model.substrategies import user_strategies

logger = custom_logger.get_logger(__name__, log_file='chem_logs.log')

class PkStrategyGenerator:
    """
    This class defines the equivalent MØD strategy for a given PKS-Domainstring. This MØD strategy is build
    iteratively depending on the contained Domains and substrates in the provided PKS. The resulting polyketide
    strategy should enable MØD to create a derivation graph, which contains the polyketide molecule that would
    otherwise be created during the polyketide synthase (PKS).
    """

    def __init__(self, assembly_line, molecules, reactions, substrategies):
        self.assembly_line = assembly_line
        self.molecules = molecules
        self.pk_reactions = reactions['fwd']['main_chain']
        self.pk_substrategies= substrategies['fwd']

        self.indigo = Indigo()

        # Counters and boundaries
        self.idx = 0
        self.assembly_items = len(assembly_line)

        # List containing performed reactions
        self.performed_reactions = []

        # Shortcuts
        self.remove_nonsubset_except_basics = _mod.filterUniverse(self._is_in_subset_or_in_basics)
        self.add_basics = _mod.addUniverse(
            molecules['basic_units']['Proton'],
            molecules['basic_units']['Hydride'],
            molecules['basic_units']['Water'],
        )


    def _is_in_subset_or_in_basics(self, graph, state, is_first_call):
        """
        Helper function to check if a graph (which represents a molecule) is in a subset 
        or in the list of basic molecules.
        Args:
            graph: The graph to check.
            state: The current state of the synthesis.
            is_first_call (bool): Indicates if it's the first call.
        Returns:
            bool: True if the graph is in the subset or in basic molecules, False otherwise.
        """
        return graph in state.subset or graph in self.molecules['basic_units'].values()


    # Functions to let the strategy grow
    def _find_molecule(self, candidates, default) -> _mod.Graph:
        """
        This function searches for a molecule within a list of candidates until it encounters
        'ACP' (Acyl Carrier Protein) in the PKS assembly line. It returns the found molecule 
        or a default molecule if no molecule is found prior to 'ACP'.
        Args:
            candidates (dict): A dictionary of candidate molecules to search in.
            default (str): The default molecule to return if no molecule is found before 'ACP'.
        Returns:
            _mod.Graph: The found molecule or the default molecule.
        """
        requested_mol = self.assembly_line[self.idx+1]
        if requested_mol in candidates:
            return candidates[requested_mol]
        else:
            message = f"Warning: Did not encounter molecule prior to 'ACP'. Using '{default}' instead."
            logger.warning(message)
            return candidates[default]


    def _find_release_substrategies(self, candidates, default):
        """
        This function iterates through the PKS assembly line from the current index
        until the end, searching for release substrategies within a list of candidate
        substrategies. It returns a list of found substrategies or a list containing the
        default substrategy if no specific substrategies are found.
        Args:
            candidates (dict): A dictionary of candidate release substrategies to search in.
            default (str): The default release substrategy to use if none are found.
        Returns:
            list: A list of found release substrategies or a list containing the default
                release substrategy.
        """
        found = []

        requested_release_type = self.assembly_line[self.idx+1]
        if requested_release_type in candidates:
            found.append(candidates[requested_release_type])

        if not found:
            # If no specific release substrategy is found, use the specified default one.
            found.append(candidates[default])

        self.performed_reactions.append('Release reactions: ' + requested_release_type)
        return found


    def _pred_ensure_min_ring_size_6(self, derivation):
        """
        Detect rings of size smaller than 6. If present return False to remove molecule.
        """
        for graph in derivation.right:
            if graph != self.molecules['basic_units']['Anchor_of_growing_chain']:
                relevant_product = graph
                break
        # Detect presence of ring with size < n with Indigo
        # Time ~50µs, 5x faster and more robust than with RDKit
        try:
            m = self.indigo.loadMolecule(relevant_product.smiles)
            if list(m.iterateRings(2, 5)):
                return False
        except Exception:
            pass
        return True


    def _pred_prevent_bridged_rings(self, derivation):
        """
        This function checks if there are bridged rings in a chemical
        derivation, specifically in the right side of the derivation. Bridged rings are
        rings where one or more atoms are shared between two or more non-fused rings. If
        bridged rings are detected, the function returns False to prevent the molecule
        from being included in the synthesis.
        Args:
            derivation: The chemical derivation to be checked for bridged rings.
        Returns:
            bool: True if no bridged rings are found, indicating that the molecule is
                allowed; False if bridged rings are detected, indicating that the
                molecule should be prevented.
        """
        for graph in derivation.right:
            if graph not in self.molecules['basic_units'].values():
                relevant_product = graph
                break
        # TODO: aromaticity seems to be wrongly encoded, rdkit fails with valence error
        ##smi = _ind.smi_to_smi_aromatic(relevant_product.smiles) # no Indigo used ATM, ignore aromaticity
        m = Chem.MolFromSmiles(relevant_product.smiles)
        num = Descriptors.rdMolDescriptors.CalcNumBridgeheadAtoms(m)
        if num > 0:
            return False
        return True


    def _pred_ensure_realism(self, derivation):
        """
        Unused function.
        """
        for graph in derivation.right:
            if graph not in self.molecules['basics'].values():
                relevant_product = graph
                break
        try:
            molecule = Chem.MolFromSmiles(relevant_product.smiles)
            Chem.rdDistGeom.EmbedMolecule(molecule, maxAttempts=1)
            return True
        except Exception:
            return False    


    def add_pk_molecule(self, strategy):
        """
        Adds a strategy containing a molecule, which is either a starter or extender molecule,
        and adds it to the Subset of available molecules.
        For more see: https://jakobandersen.github.io/mod/dgStrat/index.html#add-subset
        """
        # If we are at the very beginning, we have to pick a molecule from the starter set
        if self.idx == 0:
            molecule = self._find_molecule(self.molecules['starter_units'], 'Acetyl')

            self.performed_reactions.append('Add ' + molecule.name)
            return strategy >> _mod.addSubset(molecule)
        else:
            # Note: This AT domain is used to extend the growing PK chain,
            # in which case, we must also perform a condensation reaction (which corresponds to a KS Domain)
            molecule = self._find_molecule(self.molecules['extender_units'], 'Malonyl')

            self.performed_reactions.append('Add ' + molecule.name)
            return self.run_pk_condensation(strategy >> _mod.addSubset(molecule))


    def reduce_universe(self, strategy):
        """
        Adds a strategy to reduce the Universe (chemical space) by removing all non Subsets 
        except for when they are basic molecules (like Water).
        For more see: https://jakobandersen.github.io/mod/dgStrat/index.html#filter-universe
        """
        self.performed_reactions.append('Reduce chem. space')
        return strategy >> self.remove_nonsubset_except_basics


    def run_pk_condensation(self, strategy):
        """
        Adds a strategy to elongate the growing PK chain by including a claisen condensation reaction.
        Corresponds to a 'KS' in the PKS Domain chain.
        """
        self.performed_reactions.append('Claisen condensation')
        return strategy >> self.pk_reactions['claisen_condensation']


    def run_pk_carbonyl_reductase(self, strategy):
        """
        Adds a strategy to allow for the creation of a hydroxyl group (OH) in the PK chain.
        Corresponds to a 'KR' in the PKS Domain chain.
        """
        self.performed_reactions.append('Beta carbonyl reduction')
        return strategy >> self.pk_reactions['beta_carbonyl_reduction']


    def run_pk_dehydration(self, strategy):
        """
        Adds a strategy for removal of Water molecules in the PK chain.
        Corresponds to a 'DH' in the PKS Domain chain.
        """
        self.performed_reactions.append('Beta hydroxyl dehydration')
        return strategy >> self.pk_reactions['beta_hydroxyl_dehydration']


    def run_pk_enoyl_reductase(self, strategy):
        """
        Adds a strategy to turn a C=C double bond into a C-C single bond. (Unsaturated
        enoyl group becomes a saturated alkyl group)
        Corresponds to a 'ER' in the PKS Domain chain.
        """
        self.performed_reactions.append('Enoyl reduction')
        return strategy >> self.pk_reactions['enoyl_reduction']


    def run_release(self, strategy):
        """
        Generate a release strategy for the polyketide.

        This function constructs a release strategy for the polyketide based on various
        release mechanisms such as hydrolysis, macrolactonization, macrolactamization, and
        claisen condensation. The choice of release mechanism depends on the available
        substrategies and the specific PKS Domain string.

        Args:
            strategy: The current chemical synthesis strategy.

        Returns:
            strategy: The updated chemical synthesis strategy with the chosen release
                    mechanism added to it.
        """
        _parallel = _mod.DGStrat.makeParallel
        _rule = _mod.DGStrat.makeRule
        hydrolysis = self.pk_substrategies['release_hydroxyl_group_water']
        macrolactonization = _mod.rightPredicate[self._pred_ensure_min_ring_size_6](
            self.pk_substrategies['release_hydroxyl_group_intramolecular'])
        macrolactamization = _mod.rightPredicate[self._pred_ensure_min_ring_size_6](
            self.pk_substrategies['release_amino_group_intramolecular'])
        claisen_condensation = self.pk_substrategies['release_claisen']
        #_parallel([
        #    _rule(self.pk_reactions['release']['claisen_condensation_rs6_nc']),
        #    _rule(self.pk_reactions['release']['claisen_condensation_rs6_dc']),
        #]) # TODO: non-consuming aromatization >> _mod.revive(pk_substr['aromatization'])
        any_release = _parallel([
            hydrolysis,
            macrolactonization,
            macrolactamization,
            claisen_condensation,
        ])

        # TODO: constrainShortestPath does not work because the atoms are not connected on the lhs
        # TODO: parallel all (or better write down multiple ones, e.g. to allow claisen+hydrolysis followed by aldol?)
        candidates = {
            'Any': any_release,
            'Claisen': claisen_condensation,
            'Hydrolysis': hydrolysis,
            'Macrolactamization': macrolactamization,
            'Macrolactonization': macrolactonization,
        }

        found = self._find_release_substrategies(candidates, 'Hydrolysis')
        if len(found) == 1:
            release_strategy = found[0]
        else:
            release_strategy = _parallel(found)
        
        return strategy >> release_strategy


    def run_cyclization(self, strategy):
        """
        Generate a cyclization strategy for the polyketide synthesis.

        This function constructs a cyclization strategy for the PKS based on the
        presence of specific indicators in the PKS assembly line. This part of the strategy
        focuses on creating cyclic structures in the synthesized molecule. The function evaluates the
        text indicators to determine the type and number of cyclization events to perform.

        Args:
            strategy: The current chemical synthesis strategy.

        Returns:
            strategy: The updated chemical synthesis strategy with the cyclization steps included.
        """
        cyc_and_aro = _mod.rightPredicate[self._pred_prevent_bridged_rings](  #[_pred_ensure_realism](
            _mod.revive(self.pk_substrategies['cyclization'])
            >> _mod.revive(self.pk_substrategies['aromatization'])
        )

        try:
            requested_cyclization_number = self.assembly_line[self.idx+1]
        except IndexError:
            message = "Error parsing number argument for cyclization => ignoring cyclization all together"
            logger.warning(message)
            pass

        if requested_cyclization_number == '*':
            strategy = strategy >> _mod.repeat(cyc_and_aro)
        elif requested_cyclization_number == '0':
            pass
        else:
            try:
                number = int(requested_cyclization_number)
                strategy = strategy >> _mod.repeat[number](cyc_and_aro)
            except ValueError:
                message = "Error parsing number argument for cyclization => ignoring cyclization all together"
                logger.warning(message)
                pass
        self.performed_reactions.append('Number of cyclizations: ' + requested_cyclization_number)
        return strategy


    def run_hydroxyl_tailoring(self, strategy):
        # TODO (Robert): to be implemented?
        logger.warning("Substrategy for hydroxyl tailoring called, but not yet implemented.")
        return strategy


    def convert_assembly_line_to_strategy(self):
        """
        Convert a phenotype assembly instruction into a chemical synthesis strategy by applying the corresponding substrategy
        for each used Domain and substrate of the PKS (in the correct chronological order).

        Returns:
            DGStrat: A chemical synthesis strategy.
            List: A list containing the performed chemical reactions (in chronological order)
        """
        # Switch table for efficient lookup of text-to-function association
        # - https://jaxenter.com/implement-switch-case-statement-python-138315.html
        strategy_mapping = {
            'AT': self.add_pk_molecule, #1 # either only adds starter mol OR adds extender mol + condensation reaction
            'ACP': self.reduce_universe, #2
            'KR': self.run_pk_carbonyl_reductase, #5
            'DH': self.run_pk_dehydration, #6
            'ER': self.run_pk_enoyl_reductase, #7
            'TE': self.run_release, #8
            'Cyclization': self.run_cyclization, #9 ## comment this out during debugging to increase speed
        }

        # Load extended rule_to_func_mapping with user-defined strategies
        extended_strategy_mapping = user_strategies.extended_strategy_mapping

        # Start the strategy by adding some basic molecules
        strategy = self.add_basics

        # Iterate over the instructions of how to build the phenotype in order to choose the appropriate (sub-)strategies
        while self.idx < self.assembly_items:
            assembly_item = self.assembly_line[self.idx]

            # Check if we have a matching strategy for this part of the PKS assembly line
            # Note: We first check the user defined functions, which gives them the possibility to potentially override an existing mapping
            if assembly_item in extended_strategy_mapping:
                # We found a user defined mapping -> extend the Strategy accordingly (and provide context from this class as well)
                strategy = extended_strategy_mapping[assembly_item](strategy, self)
            elif assembly_item in strategy_mapping:
                # We found a predefined mapping -> extend the Strategy accordingly
                strategy = strategy_mapping[assembly_item](strategy)
            
            self.idx += 1 # Increment index and continue looking for matches until end of PKS assembly line
            # TODO (Robert): alpha modifications

        # TODO: Temporarilly added post_pks modifications (but doesn't seem to affect generated molecules)
        strategy >> self.pk_substrategies['post_pks']
 
        return strategy, self.performed_reactions