from typing import List

import config.custom_logger as custom_logger
from config.pkevo_config import SearchConfig as search_config
from grammatical_evolution.model._individual import _Individual

logger = custom_logger.get_logger(__name__, "ge_logs.log")

class FitnessFunctions:
    class FunctionTemplate:
        """
        This servers as a boilerplate template for writing custom fitness functions. The way it works is that you
        write your fitness function as a class (give it an arbitrary, yet meaningful name) inside this module and
        define the classattributes and write your fitness function inside 'objective_function' (the name needs to 
        remain unchanged!). You can of course add additional functions to your class, but the application will only
        ever call the 'objective_function' function.
        """
        # (OPTIONAL) add import statements here or inside your functions (for performance reasons, don't write your imports 
        # at the top of the file, so that they are only ever imported if your fitness function is actually going to be used)

        # Class variables (no need for instance variables, since these variables will get passed onto them)
        coevolution: bool # if you want to evaluate the whole population at once, set this to 'True', use 'False' otherwise
        default_fitness = float # define what the default fitness of an unevaluated individual should be
        optimization_type = str # use either 'min' or 'max' depending on what you are trying to achieve
        is_chemistry_necessary = bool # if 'True', MØD will be used in the GE process in order to generate molecules, which can be your basis for evaluation
        target = str # you need to define your target (towards which the GE alg. will try to evolve) inside the SearchConfig class (/config/pkevo_config.py)
        
        # Instance method(s)
        def objective_function(self, individual_s, cached_dict) -> float | None:
            """
            This is the heart and brain of the whole class. What you do here is completely up to you. The only condition 
            is that you return a float value (which will be set as the fitness score of the individual) if you're 
            evaluating Individuals separately (when `coevolution` is set to False), otherwise when evaluating whole 
            populations at once, you do not return anything, but instead set the fitness value for each individual 
            directly inside this function (since Individuals are passed as references, those changes will remain 
            even outside the function). You may even choose to cause "side-effects" by editing the Individual further 
            (e.g. giving its molecules individual fitness values).

            Args:
                individual_s (Individual | List[_Individual]): Depending on value of`coevolution`, either single Individual 
                    or list of Individuals
                cached_dict (dictionary): you may use this to increase performance of your function by storing and looking 
                    up already calculated results
            Returns:
                float | None: Return the fitness value when evaluating single Individuals, otherwise don't return anything

            """
            # Note: Technically, it wouldn't be necessary to return the calculated fitness value, due to the Individuals being 
            # passed as references. However, I consider it crucial to point out that the goal of the objective function is to 
            # assign a fitness value to an individual, which in turn is one of the (if not the) most important concepts behind 
            # evolutionary algorithms. Hence, the explicit return of the fitness value.

            # Note 2: However, due to the `coevolution` already kind of breaking my approach pointed out above, this behaviour 
            # of explicit returns might be removed in the future.

            # Note 3: Sorting the Individual.molecules list is not necessary. It will be sorted by the GE algorithm afterwards 
            # from best to worst based on the optimization type of your fitness function
            pass


    class StringDistance:
        coevolution = False
        default_fitness = 100000.0
        optimization_type = 'min'
        is_chemistry_necessary = False
        target = search_config.TARGET_SEQUENCE

        # custom constructor only necessary to setup a Levensthein implementation
        def __init__(self):
            try:
                from Levenshtein import distance as lev
                self.lev = lev  # Assign the imported function to an instance attribute
            except ImportError:
                self.lev = None  # Define lev as None if the import fails


        def _levenshtein_distance(self, s1, s2):
            """
            Compute the Levenshtein distance between two strings, which is the minimum number of single-character edit operations
            (insertions, deletions, or substitutions) required to transform one string into another.

            Parameters:
                s1 (str): The first string.
                s2 (str): The second string.

            Returns:
                int: The Levenshtein distance between the two strings.
            """
            if self.lev is not None:
                # use the 3rd party library 'Levenshtein'
                return self.lev(s1, s2)
            else:
                # use previous implementation by Robert Haas
                if len(s1) < len(s2):
                    return self._levenshtein_distance(s2, s1)

                if len(s2) == 0:
                    return len(s1)

                previous_row = range(len(s2) + 1)
                for i, c1 in enumerate(s1):
                    current_row = [i + 1]
                    for j, c2 in enumerate(s2):
                        insertions = previous_row[j + 1] + 1
                        deletions = current_row[j] + 1
                        substitutions = previous_row[j] + (c1 != c2)
                        current_row.append(min(insertions, deletions, substitutions))
                    previous_row = current_row

                return previous_row[-1]


        def objective_function(self, individual: _Individual, remembered_distances: dict):
            """
            Calculate the fitness of the generated phenotype (= PKS Domain String) based on its Levenshtein distance from the target sequence.

            Parameters:
                generated_sentence (str): The sentence generated by the GE algorithm.
                remembered_distances (dict): A dictionary to store already calculated distances for caching.

            Returns:
                int: The fitness of the generated sentence based on the Levenshtein distance.
            """
            target_sequence = self.target
            generated_sentence = individual.phenotype

            if generated_sentence in remembered_distances:
                distance = remembered_distances[generated_sentence]
            else:
                distance = self._levenshtein_distance(generated_sentence, target_sequence)
                remembered_distances[generated_sentence] = distance

            return distance


    class FingerprintSimilarity:

        coevolution = False
        default_fitness = 0.0
        optimization_type = 'max'
        is_chemistry_necessary = True
        target = search_config.TARGET_SMILES


        def objective_function(self, individual: _Individual, remembered_similarities: dict):
            """
            Calculating molecular similarity using fingerprints (in this case using Tanimoto similarity) as shown 
            here: https://www.rdkit.org/docs/GettingStartedInPython.html#fingerprinting-and-molecular-similarity
            """

            from rdkit import Chem
            from rdkit.Chem import AllChem, DataStructs

            target_smiles = self.target

            # set defaults
            max_similarity = 0.0

            if len(individual.molecules) == 0:
                # if applied strategy was not able to generate a molecule, we have no SMILES string => no similarity
                return max_similarity

            # Convert target smiles to a molecule
            target_molecule = Chem.MolFromSmiles(target_smiles)

            # Calculate (Morgan) fingerprint for the target molecule
            radius = 2  # You can adjust the radius based on your preference
            #fpgen = AllChem.GetRDKitFPGenerator()
            #target_fp = fpgen.GetFingerprint(target_molecule)
            target_fp = AllChem.GetMorganFingerprintAsBitVect(target_molecule, radius)

            # Calculate similarities for each generated molecule against target molecule
            for mol_obj in individual.molecules:

                # perform lookup to see if this molecule was generated (and therefore evaluated) before
                if mol_obj.smiles_string in remembered_similarities:
                    mol_obj.fitness = remembered_similarities[mol_obj.smiles_string]
                else:
                    # molecule is unknown -> run fingerprint similarity alg from RDKit and store similarity result in our Molecule object
                    generated_molecule = Chem.MolFromSmiles(mol_obj.smiles_string)

                    # Skip if generated molecule cannot be parsed, however store this 'failed' molecule in our dictionary, 
                    # so that we don't keep repeating this 'failure' over and over again
                    if generated_molecule is None:
                        remembered_similarities[mol_obj.smiles_string] = 0.0
                        mol_obj.fitness = 0.0
                        continue 

                    # generated_fp = fpgen.GetFingerprint(generated_molecule)
                    generated_fp = AllChem.GetMorganFingerprintAsBitVect(generated_molecule, radius)
                    # possibly good alternative similarity algorithms for Tanimoto are Sokal and Kulczynski
                    similarity = DataStructs.TanimotoSimilarity(generated_fp, target_fp)

                    # store the similarity result as the fitness of the Molecule object of the Individual
                    mol_obj.fitness = similarity
                    # also store the result in the overall lookup dictionary
                    remembered_similarities[mol_obj.smiles_string] = similarity

                if mol_obj.fitness > max_similarity:
                    max_similarity = mol_obj.fitness      

            return max_similarity


    class PharmacophoreScore:

        coevolution = True
        default_fitness = 0.0
        optimization_type = 'max'
        is_chemistry_necessary = True
        target = search_config.TARGET_PHARMACOPHORE_FILE


        def objective_function(self, individuals: List[_Individual], remembered_pharma_scores: dict):

            from external_tools.ligand_scout_singleton import ligand_scout_utils

            pharmacophore_file = self.target
            score_list = []

            # get all unique SMILES strings from individuals
            new_molecules_dict = dict()
            for ind in individuals:
                if len(ind.molecules) > 0:
                    for mol in ind.molecules:
                        # check if we have seen this molecule before
                        if mol.smiles_string not in remembered_pharma_scores:
                            new_molecules_dict[mol.smiles_string] = 0.0 # initialize fitness of molecule with 0
                            # Initialize this new molecule in the cache dictionary
                            remembered_pharma_scores[mol.smiles_string] = 0.0
                            
            # extract only the SMILES representation of the new molecules into a list
            new_smiles_list = list(new_molecules_dict.keys())

            if len(new_smiles_list) == 0:
                # there are no unknown molecules which need to be evaluated by LigandScout
                logger.warning(f"Generation didn't generate a single new molecule.")
            else:
                #print(f"{len(new_smiles_list)=}")
                score_list = ligand_scout_utils.get_pharmacophore_score_list(new_smiles_list, pharmacophore_file, min_config=True)

            # Use zip to combine the LigandScout score results of new molecules back to the dictionaries
            for smiles, score in zip(new_molecules_dict.keys(), score_list):
                new_molecules_dict[smiles] = score
                remembered_pharma_scores[smiles] = score

            for ind in individuals:
                if len(ind.molecules) > 0:
                    overall_best_fitness = 0.0
                    for mol in ind.molecules:
                        # Note: don't use the new_molecules dict since an Individual can have 'old' molecules as well!
                        mol.fitness = remembered_pharma_scores[mol.smiles_string]
                        if mol.fitness > overall_best_fitness:
                            overall_best_fitness = mol.fitness
                    ind.fitness = overall_best_fitness # Individual.fitness must be set here, since we skipped the default Individual.evaluate() function!


