# Based on Robert Haas' _Grammar class from "panakeias_garden" package (/panakeias_garden/pre_release_demos/metaheuristics.py)
# and PonyGE's Grammar class (https://github.com/jmmcd/ponyge/tree/master) 

import io
import re
import time

from config.file_storage_singleton import file_storage


class _Grammar():
    """
    A GE Context-Free Grammar.

    This class represents a context-free grammar used in the Grammatical Evolution algorithm. It parses
    the BNF grammar text and provides methods to generate phenotypes from genomes using the grammar rules.

    Attributes:
        NT (str): Non-terminal symbol marker ("<NT>"). Used to distinguish non-terminals in the grammar.
        T (str): Terminal symbol marker ("<T>"). Used to distinguish terminals in the grammar.
        rules (dict): A dictionary mapping non-terminal symbols to their corresponding production choices.
        non_terminals (set): A set containing all non-terminal symbols present in the grammar.
        terminals (set): A set containing all terminal symbols present in the grammar.
        start_rule (tuple): A tuple representing the start rule of the grammar in the form (start_symbol, type).
                            'start_symbol' is the non-terminal representing the starting rule, and 'type' is 'NT'.

    Note:
        The _Grammar class is designed to be used internally by the GrammaticalEvolution algorithm. It provides
        methods for generating phenotypes and should not be accessed directly outside the algorithm implementation.
    """

    NT = "NT" # Non-Terminal
    T = "T" # Terminal


    def __init__(self, grammar_text=None, grammar_file_path=None, verbose=False, store_results=False):
        """
        Constructor for the _Grammar class.

        Args:
            grammar_text (str, optional): The BNF grammar text. Defaults to None.
            grammar_file_path (str, optional): The path to a file containing the BNF grammar text. Defaults to None.
            verbose (bool, optional): Whether to print verbose output. Defaults to False.
            store_results (bool, optional): Whether to store results as files. Defaults to False.
        """
        print("\n=======================")
        print("Setting up grammar...")
        print("=======================\n")

        start_time = time.perf_counter()

        self.rules = {}
        self.non_terminals, self.terminals = set(), set()
        self.start_rule = None

        if not grammar_text:
            if not grammar_file_path:
                raise ValueError("Either grammar_text or grammar_file_path must be provided.")
            print(f"Using grammar defined in {grammar_file_path}")
            grammar_text = self.read_grammar_file(grammar_file_path)

        self._parse_grammar_text(grammar_text)

        end_time = time.perf_counter()
        elapsed_time = end_time - start_time

        if verbose:
            print(self)

        if store_results:
            file_storage.save_text_file("grammar.txt", self.__str__())


        print(f"\nLoading grammar took (seconds): {elapsed_time:.3f}")
        print()
        print("=======================")
        print("Grammar setup complete.")
        print("=======================")
        print()


    def __str__(self):
        """
        Returns a string representation of the _Grammar object, showing non-terminals, terminals, and rules.

        Returns:
            str: A string representation of the grammar.
        """
        header = f"Grammar with following {len(self.non_terminals)} non-terminals, {len(self.terminals)} terminals and {len(self.rules)} rules."
        non_terminal_info = "\n\nNon-terminals:\n" + ", ".join(self.non_terminals)
        terminal_info = "\n\nTerminals:\n" + ", ".join(self.terminals)
        rule_info = self.format_rules()
        
        return header + non_terminal_info + terminal_info + rule_info


    def format_rules(self):
        """
        Formats the grammar rules for display.

        Returns:
            str: A string representation of the rules.
        """
        rules = []
        for index, (k, v) in enumerate(self.rules.items(), start=1):
            rules.append(f"  {index}. Rule: {k} -> {v}")
        rule_info = "\n\nRules:\n" + "\n".join(rules) + "\n"
        return rule_info


    def read_grammar_file(self, grammar_file_path):
        """
        Read a grammar file from the given path.

        Args:
            grammar_file_path (str): The path to the grammar file.

        Returns:
            str: The contents of the grammar file as a string.
        """
        try:
            # Read the grammar file
            with open(grammar_file_path, 'r') as file:
                grammar_text_bnf = file.read()
            return grammar_text_bnf
        except FileNotFoundError:
            print(f"Grammar file '{grammar_file_path}' not found.")
            exit(1)


    def _parse_grammar_text(self, grammar_text):
        """
        Parse the BNF grammar text.

        Args:
            grammar_text (str): The BNF grammar text.
        """
        non_terminal_pattern = "(<.+?>)"   # regex <.+?> Non greedy match of anything between brackets => condition: NTs must be enclosed by <>
        rule_separator = "::="
        production_separator = "|"

        # Read the grammar file
        for line in io.StringIO(grammar_text):
            if not line.startswith("#") and line.strip() != "":
                # Split rules. Everything must be on one line
                if line.find(rule_separator):
                    lhs, productions = line.split(rule_separator)
                    lhs = lhs.strip()
                    if not re.search(non_terminal_pattern, lhs):
                        raise ValueError("lhs is not a NT in the form: (<.+?>)", lhs)
                    self.non_terminals.add(lhs)
                    if self.start_rule == None:
                        self.start_rule = (lhs, self.NT)
                    # Find terminals
                    tmp_productions = []
                    for production in [production.strip()
                                       for production in
                                       productions.split(production_separator)]:
                        tmp_production = []
                        if not re.search(non_terminal_pattern, production):
                            self.terminals.add(production)
                            tmp_production.append((production, self.T))
                        else:
                            # Match non terminal or terminal pattern
                            for value in re.findall("<.+?>|[^<>]*", production):
                                if value != '':
                                    if not re.search(non_terminal_pattern,
                                                     value):
                                        symbol = (value, self.T)
                                        self.terminals.add(value)
                                    else:
                                        symbol = (value, self.NT)
                                    tmp_production.append(symbol)
                        tmp_productions.append(tmp_production)
                    # Create a rule
                    if not lhs in self.rules:
                        self.rules[lhs] = tmp_productions
                    else:
                        raise ValueError("lhs should be unique", lhs)
                else:
                    raise ValueError("Each rule must be on one line")
    

    def generate_phenotype(self, genome, max_wraps=2):
        """
        Generate a phenotype from the given genome using the grammar rules.

        Parameters:
            genome (list): The genome (sequence of integers) to be used for generating the phenotype.
            max_wraps (int, optional): The maximum number of times the genome can wrap during the generation process.
                                   Wrapping occurs when the genome reaches its end and needs to restart
                                   from the beginning for further symbol expansion.

        Returns:
            tuple: A tuple containing the generated phenotype, the number of used codons in the genome, and the PKS assembly items.
                The phenotype is a string representation of the generated output.
                The used_codons is an integer representing the number of codons consumed.
                The assembly_items is a list of domains and molecules that comprise the PKS.
        """
        used_codons = 0 # Tracks how much of the input genome has been used in the mapping
        wraps = 0 # Tracks the number of wraps (cycles) performed
        phenotype = [] # This stores the resulting phenotype
        used_symbols = []  # Store the used Non-terminals and Terminals (in occuring order) (This is more informational than `assembly_items`)
        assembly_items = [] # Store the domains and molecules that comprise the PKS (in occuring order)
        production_choices = [] # Stores the available production choices for expansion

        # start the list of unexpanded symbols with the start rule `<S>`
        unexpanded_symbols = [self.start_rule]

        # As long as we haven't exceeded our wrapping limit and we still have unexpanded symbols left, continue with mapping the genome
        while (wraps <= max_wraps) and (len(unexpanded_symbols) > 0):
            
            # Wrap around the genome if we reached its end and we are still left with production choices (which indicates we have NTs remaining)
            # Spoken in simple terms: Go back to start of genome and continue the mapping process if some NTs remained unmapped
            if used_codons % len(genome) == 0 and used_codons > 0 and len(production_choices) > 1:
                wraps += 1

            # Expand a production by taking the first symbol from unexpanded_symbols
            current_symbol = unexpanded_symbols.pop(0)

            # If the current symbol is a terminal, add it to the phenotype and store the used symbol (= the RHS of rule applied)
            if current_symbol[1] == self.T:
                production_choices = [] # Terminals offer no production choices. Important to reset to empty for when wrapping condition is checked
                phenotype.append(current_symbol[0])
                if current_symbol[0] != '-' and current_symbol[0] != '--' and current_symbol[0] != ')':
                    assembly_items.append(current_symbol[0].strip("("))
                used_symbols.append(current_symbol[0])  # Store the value of the used rule
            else:
                # The current symbol is a non-terminal -> get the production choices for it
                production_choices = self.rules[current_symbol[0]]
                used_symbols.append(current_symbol[0])  # Store the name of the used rule

                # Select a production choice based on current position on genome (`% len(genome)` needed due to possible wrapping)
                # and number of current choices `len(production_choices)`
                current_production_position = genome[used_codons % len(genome)] % len(production_choices)

                # Increment used_codons if there was more than 1 choice
                if len(production_choices) > 1:
                    used_codons += 1

                # Add the symbols of the selected production choice at the start of unexpanded_symbols (Depth-first approach)
                unexpanded_symbols = production_choices[current_production_position] + unexpanded_symbols

        # If there are unexpanded symbols left even after wrapping x times, the mapping is not complete -> Individual has no phenotype and is invalid
        if len(unexpanded_symbols) > 0:
            return (None, used_codons, None)
        phenotype = "".join(phenotype)
        return (phenotype, used_codons, assembly_items)


    def deconstruct_string(self, input_string):
        """
        Extract valid terminal symbols from an input string and return the PKS assembly line.

        This is just a helper function not used in the general PKevo workflow. Its purpose is to
        quickly translate a PKS Domain string into an assembly line, which can be used by MØD.

        Caution: This function does not guarantee a PKS to be valid in terms of the used grammar rules.
        It only takes into account the terminals from that grammar and nothing else. So it is possible to
        get PKSs which would otherwise not be allowed by the selected grammar. Use this function at 
        your on descretion!

        Args:
            input_string (str): The input string containing terminal symbols enclosed in parentheses.

        Returns:
            list: A list of cleaned-up terminal symbols found in the input string.
        """
        
        # Sort the terminals list by length in descending order to prioritize longer matches
        sorted_terminals = sorted(self.terminals, key=len, reverse=True)
        # Turn terminals into a search pattern for regex, but exclude certain terminals like dashes
        pattern = f'({"|".join(re.escape(sym) for sym in sorted_terminals if sym not in {"-", "--", ")"})})'
        # Find matches using regex and the defined patterns, but exclude empty strings
        matches = [match for match in re.findall(pattern, input_string) if match != ""]
        # Final cleanup of the results
        cleaned_matches = [match.strip("()") for match in matches]
        return cleaned_matches