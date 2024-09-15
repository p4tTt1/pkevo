import os
import time

from external_tools.mod_importer import mod
import config.custom_logger as custom_logger

from .substrategy_generator import SubstrategyGenerator
from .pk_strategy_generator import PkStrategyGenerator

logger = custom_logger.get_logger(__name__, log_file='chem_logs.log')


class Chemistry:
    """
    The Chemistry class represents the chemical environment for the polyketides and provides
    methods for loading molecules, reactions, and generating polyketide synthesis strategies.

    Attributes:
        molecules (dict): A dictionary containing loaded molecular structures for starter
            units, extender units, and basic units.
        reactions (dict): A dictionary containing loaded chemical reactions categorized as
            forward ('fwd') and reverse ('rev') strategies.
        release_reactions (list): Subset of `reactions` containing only those reactions involved
            in releasing the polyketides
        substrategies (dict): A dictionary of generated substrategies for polyketide synthesis.
    """

    def __init__(self, verbose = True):
        """
        Initializes the Chemistry object, loading molecules, reactions, and generating substrategies.
        """
        # Get the absolute path of the directory containing the current script file
        self.curr_path = os.path.dirname(os.path.abspath(__file__))

        start_time = time.perf_counter()
        self.molecules = {
            'starter_units': self._load_molecules("starter_units"),
            'extender_units': self._load_molecules("extender_units"),
            'basic_units': self._load_molecules("basics"),
        }
        end_time = time.perf_counter()
        elapsed_time_mols = end_time - start_time
        
        start_time = time.perf_counter()
        self.reactions = self._load_reactions()
        self.final_cascade_reactions = self._extract_final_cascade_reactions()
        end_time = time.perf_counter()
        elapsed_time_reacts = end_time - start_time
        
        start_time = time.perf_counter()
        # Note: substrategies generation (just like applied_strategy) is moved to separate script file to keep things tidy
        # self.substrategies = subs.generate_substrategies(self.molecules, self.reactions) # old script should no longer be used
        sub_generator = SubstrategyGenerator(self.molecules, self.reactions) # Instead use the new SubstrategyGenerator class
        self.substrategies = sub_generator.generate_substrategies()
        end_time = time.perf_counter()
        elapsed_time_strats = end_time - start_time

        if verbose:
            print()
            print("=======================")
            print("Setting up chemistry...")
            print("=======================")
            print()
            print(f"Loading molecules took (seconds): {elapsed_time_mols:.3f}")
            print(f"Loading reactions took (seconds): {elapsed_time_reacts:.3f}")
            print(f"Loading substrategies took (seconds): {elapsed_time_strats:.3f}")
            print()
            print("=======================")
            print("Chemistry setup complete.")
            print("=======================")
            print()


    def __reduce__(self):
        # Note: This is necessary for pickling the Chemistry class, which is needed for multiprocessing.
        # We return a tuple with a callable (constructor) and its arguments. In this case, 
        # we're using the __init__() method to reconstruct the object.
        return (self.__class__, (False,))


    def _get_file_content(self, filepath):
        """
        Helper method to read and return the content of a file given its filepath.
        
        Returns:
            str: The content of the file as a string
        """
        with open(filepath) as file_handle:
            content = file_handle.read()
        return content


    def _load_molecules(self, target_dir):
        """
        Loads molecular structures from a specified directory and returns them as a dictionary.

        Returns:
            dict: The dictionary containing the molecules as `mod.Graph` objects
        """
        molecule_dict = {}
        units_folder = os.path.join(self.curr_path, "chemistry_model", "molecules", target_dir)
        for full_filename in os.listdir(units_folder):
            if full_filename.endswith(".smi"):
                filepath = os.path.join(units_folder, full_filename)
                content = self._get_file_content(filepath)
                content = content.strip(' \t\n\r')
                molecule_name = full_filename[:-4]  # Remove the last 4 characters (".smi") to exclude file extension
                try:
                    graph = mod.smiles(content, name=molecule_name)
                    molecule_dict[molecule_name] = graph
                except mod.InputError:
                    message = ('The content of following SMILES file could not be converted into a '
                            'molecule, therefore it was skipped: {}'.format(filepath))
                    logger.warning(message)
            else:
                logger.warning('The file extension of following file is unexpected, therefore it was skipped: {}.'.format(full_filename))
        return molecule_dict


    def _extract_final_cascade_reactions(self):
        """ 
        This method collects all reactions that can take part in the final cascade of the
        PKS (so all reactions defined inside the 'release', 'cyclization', 'aromatization' or 'post_pks' folders). 
        
        Returns:
            list: A list containing all chemical reactions for the final PKS cascade (as `mod.Rules` objects)
        """
        release_reactions = list(self.reactions['fwd']['release'].values())
        cyclization_reactions = list(self.reactions['fwd']['cyclization'].values())
        aromatization_reactions = list(self.reactions['fwd']['aromatization'].values())
        post_pks_reactions = list(self.reactions['fwd']['post_pks'].values())
        return release_reactions + cyclization_reactions + aromatization_reactions + post_pks_reactions


    def _load_reactions(self):
        """
        Loads chemical reactions from GML files in specific directories and converts them into Rule objects using MØD. 
        Reactions are categorizes as forward and reverse reactions.
        Note: Reverse reactions (Useful for deconstructing molecules into PKS Domains) not yet implemented!

        Returns:
            dict: The dictionary containing forward and reverse chemical reactions as `mod.Rules` objects
        """
        reactions = {
            'fwd': {}, # Forward strategies are for constructing the molecules (Domain String -> Molecule)
            'rev': {}, # NOT IN USE! Reverse strategies are for deconstructing the molecules (Molecule -> Domain String)
        }
        reactions_folder = os.path.join(self.curr_path, "chemistry_model", "reactions")
        
        for reaction_dir in os.listdir(reactions_folder):
            reaction_dir_path = os.path.join(reactions_folder, reaction_dir)
            if os.path.isdir(reaction_dir_path):
                # Iterate over all files and subfolders in the specified directory
                for full_filename in os.listdir(reaction_dir_path):
                    # We are only interested in chemical reactions, which are stored in GML files
                    if full_filename.endswith('.gml'):
                        filepath = os.path.join(reaction_dir_path, full_filename)
                        content = self._get_file_content(filepath)
                        reaction_name = full_filename[:-4]  # Remove the last 4 characters (".gml") to exclude file extension
                        # Create forward reaction rules
                        try:
                            rule_fwd = mod.ruleGMLString(content)
                            if reaction_dir not in reactions['fwd']:
                                reactions['fwd'][reaction_dir] = {}
                            reactions['fwd'][reaction_dir][reaction_name] = rule_fwd                            
                        except mod.InputError:
                            logger.warning(f'Warning: Reading "{filepath}" failed. '
                                'Text could not be interpreted as graph '
                                'grammar rule and will therefore be ignored.')
                            continue
                        # Create reverse reaction rules
                        # try:
                        #     #Note: Some constraints make inversions impossible. Could this have any neg. side-effects?
                        #     rule_rev = mod.ruleGMLString(content, invert=True)
                        #     if reaction_dir not in reactions['rev']:
                        #         reactions['rev'][reaction_dir] = {}
                        #     reactions['rev'][reaction_dir][reaction_name] = rule_rev
                        # except mod.InputError:
                        #     logger.warning(f'Warning: Reading "{filepath}" failed. '
                        #         'Text could not be interpreted as inverse graph '
                        #         'grammar rule and will therefore be ignored.')
                        #     continue
        return reactions


    def create_strategy_for_ind(self, assembly_line):
        """
        Creates a synthesis strategy based on the given PKS assembly line.

        Args:
            assembly_line (List[str]): list of strings containing the assembly instructions of this PKS

        Returns:
            mod.DGStrat: The synthesis strategy for the PKS
            List: List containing the necessary chemical reactions
        """
        strat_generator = PkStrategyGenerator(assembly_line, self.molecules, self.reactions, self.substrategies)
        applied_strategy, performed_reactions = strat_generator.convert_assembly_line_to_strategy()
        return applied_strategy, performed_reactions


    def convert_strategy_to_derivation_graph(self, strategy):
        """
        Converts a synthesis strategy into a derivation graph (which represents the chemical 
        reaction network) using available molecules.
        
        Args:
            strategy (mod.DGStrat): The synthesis strategy for the PKS

        Returns:
            mod.DG: derivation graph (which represents the chemical reaction network)
        """
        known_graphs = list(self.molecules['starter_units'].values()) + \
            list(self.molecules['extender_units'].values()) + \
            list(self.molecules['basic_units'].values()) 

        settings = mod.LabelSettings( 
            mod.LabelType.Term,
            mod.LabelRelation.Unification,
            False,
            mod.LabelRelation.Unification)
        # Note from MOD changelog: dgRuleComp and DG.calc have been deprecated, and their implementation is now based on DGBuilder.execute(). Use DGBuilder.execute() directly instead.
        #derivation_graph = mod.dgRuleComp(known_graphs, strategy, settings)     
        #derivation_graph.calc()
        derivation_graph = mod.DG(graphDatabase=known_graphs, labelSettings=settings)
        derivation_graph.build().execute(strategy, verbosity=1)

        return derivation_graph


    def get_smiles_for_ind(self, assembly_line, get_reactions = False):
        """
        Extracts SMILES representations for molecules synthesized using a given PKS assembly line (list object).
        
        Args:
            assembly_line (List[str]): The synthesis strategy for the PKS
            get_reactions (boolean): If True, the list of performed reactions is returned as well (defaults to False)

        Returns:
            released_products (List[str]): List containing released molecules (in SMILES representation)
            (Optional) performed_reactions (List[str]): List containing performed chemical reactions (depends on `get_reactions`)
        """
        strategy, performed_reactions = self.create_strategy_for_ind(assembly_line) # makes (potentially costly) call to MØD
        # Note: dg is of type: <class 'mod.libpymod.DG'>
        dg = self.convert_strategy_to_derivation_graph(strategy) # makes (potentially costly) call to MØD
        #print(dg.print()) # with dg.print() I get output stored persistantly (as files in /out directory)
       
        # Note: instead of looping over dg._products (which is basically a subset of all vertices) to get the molecules as smiles I instead 
        # have to check the vertices directly and make sure that they don't have any outgoing edges (a product dg._products might 'only' be 
        # an intermediate product, which is later on used as an educt for another product. I'm not interested in such intermediate products)
        products = []
        for vertex in dg.vertices: # Note: vertex is of type: <class 'mod.libpymod.DGVertex'>
            found_molecule = False
            if vertex.inDegree > 0 and vertex.outDegree == 0:
                in_edges = vertex.inEdges
                for incoming in in_edges:
                    rules = incoming.rules
                    for rule in rules:
                        # Only if incoming edge belongs to a "final cascade" reaction, the molecule is deemed a viable product
                        if rule in self.final_cascade_reactions:
                            found_molecule = True
                            molecule = vertex.graph
                            products.append(molecule)
                            break
                    if found_molecule:
                        break
                # molecule = vertex.graph
                # products.append(molecule)
        
        # Further filter the products and ignore all basic units (Water, etc) and the single Silver atom
        released_products = [
            item.smiles for item in products
            if item not in self.molecules['basic_units'].values() and '[Ag]' not in item.smiles
        ]

        if get_reactions:
            return released_products, performed_reactions, dg
        else:
            return released_products