class EvolutionConfig:
    GRAMMAR_FILE_NAME = 'pks_simplified.bnf' # Note: changing this is potentially dangerous, since matching grammar rules to substrategies might fail!
    CODON_SIZE = 127 # defines the value range limit of what each individual codon can take
    GENOME_LENGTH = 100 # set the default length of an individual's genome
    ELITE_SIZE = 1 # Defines how many of the best individuals of old generation will be considered for moving onto the next generation
    POPULATION_SIZE = 250 # determines the number of individuals (= size of the population) in the initial generation
    GENERATION_SIZE = 250 # number of individuals that will be retained in each generation after the evolutionary process
    GENERATIONS = 50 # how many generations should be simulated? This determines the number of iterations in the evolutionary loop

    MUTATION_PROBABILITY = 0.05 # the probability of a mutation occurring on a codon
    DELETION_PROBABILITY = 0.1 # the probability of a single codon being deleted from the genome
    INSERTION_PROBABILITY = 0.1 # the probability of a single codon being inserted into the genome
    CROSSOVER_PROBABILITY = 0.9 # the probability of crossovers occurring during offspring generation

    SELECTION_STRATEGY = 'truncation_selection' # Choose the selection strategy for selecting parent individuals
    REPLACEMENT_STRATEGY = 'generational_replacement' # Choose the replacement strategy for populations
    CROSSOVER_OPERATOR = 'twopoint_crossover' # Choose the crossover operation for offspring generation
    WITHIN_CODING_REGION = True # If set True, crossover and mutations only occur in coding region of the genome
    VARIABLE_LENGTH = False # If set True, crossover points on two genomes may differ, leading to either growing or shrinking offspring genomes


class SearchConfig:
    # Options to prematurely end the evolutionary loop
    FITNESS_THRESHOLD = None # If set (as a float value), the GE algorithm will stop once the defined fitness threshold is surpassed
    BREAK_AFTER = None # If set (as an int value), the GE algorithm will stop if the overall best fitness has not improved for x GE iterations 

    # Depending on selected fitness function, the algorithm either uses TARGET_SEQUENCE, TARGET_SMILES or TARGET_PHARMACOPHORE_FILE, so you don't have to declare and update all variables
    # Target for now is Geldanamycin
    TARGET_SEQUENCE = 'AT(AHBA)-ACP--KS-AT(Methylmalonyl)-KR-DH-ER-ACP--KS-AT(Methoxymalonyl)-KR-DH-ER-ACP--KS-AT(Methylmalonyl)-KR-ACP--KS-AT(Methylmalonyl)-KR-DH-ACP--KS-AT(Methoxymalonyl)-KR-ACP--KS-AT(Malonyl)-KR-DH-ER-ACP--KS-AT(Methylmalonyl)-KR-DH-ACP--TE(Macrolactonization)'
    TARGET_SMILES = 'CC1CC(C(C(C=C(C(C(C=CC=C(C(=O)NC2=CC(=O)C(=C(C1)C2=O)OC)C)OC)OC(=O)N)C)C)O)OC'
    TARGET_PHARMACOPHORE_FILE = 'gdm_pharmacophore_relaxed.pmz'


class GeneralConfig:
    MULTIPROCESSING_CORES = 4 # Select the number of cores to be used for multiprocessing (Set it to `1` to disable multiprocessing)
    LIGANDSCOUT_PATH = "/home/wserver/wsoft/ligandscout4" # path to where PKevo can access LigandScout
    MOD_PATH = "/home/talax/xtof/local/Mod/lib64" # path to where PKevo can access MedOlDatschgerl (MOD)
    GENERATE_PLOTS = True # If set True, the app will generate plots containing statistical measures for the created generations and individuals