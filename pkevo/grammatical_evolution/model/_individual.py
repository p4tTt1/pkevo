# Based on Robert Haas' _Individual class from "panakeias_garden" package (/panakeias_garden/pre_release_demos/metaheuristics.py) 

import random
import re
from typing import List

from . import pet_names as pn
from ._pks_architecture import PKSArchitecture
from ._molecule import Molecule

from config.pkevo_config import EvolutionConfig as evo_config


class _Individual():
    """
    This class represents an individual in the population.

    Attributes:
        genome (list[int]): The genome of the individual, which is a sequence of integers.
        optimization_type (str): The type of optimization, either "max" for maximization or "min" for minimization.
        fitness (float): The fitness value of the individual, computed by the fitness function.
        phenotype (str): The phenotype generated using the grammar from the genome.
        used_codons (int): The number of codons used in the genome to generate the phenotype.
        assembly_line (object): Building instructions for the phenotype (ordered list of applied rules to generate phenotype).
        pet_name (str): A name to make it more personal (if given in init(), use it as a constraint for generating a new name).
        architecture (PKSArchitecture): Object describing the architecture of this PKS in a more human-readable way

    Note:
        The _Individual class is intended to be used as an internal representation of individuals in the
        Grammatical Evolution algorithm. It is not meant to be accessed directly outside the algorithm implementation.
    """


    def __init__(self, optimization_type, default_fitness, genome=None, length=evo_config.GENOME_LENGTH, codon_size=evo_config.CODON_SIZE, name_restriction=None):
        """
        Constructor for the _Individual class.

        Args:
            optimization_type (str): The type of optimization, either "max" for maximization or "min" for minimization.
            default_fitness (float): The default fitness value of the individual.
            genome (list[int], optional): The genome of the individual. Defaults to None.
            length (int, optional): The length of the genome. Defaults to 100.
            codon_size (int, optional): The size of codons. Defaults to 127.
            name_restriction (str, optional): Giving individuals names makes it more personal. With this its possible
                                              to restrict the creature, this individuals can become. Defaults to None.
        """

        if genome == None:
            # by default this creates a list of 100 random integers, each between 0-127
            self.genome = [random.randint(0, codon_size)
                           for _ in range(length)]
        else:
            self.genome = genome

        self.optimization_type = optimization_type
        self.fitness = default_fitness
        self.phenotype = None
        self.molecules: List[Molecule] = []  # Initialize as an empty list for storing Molecule objects
        self.used_codons = 0
        self.assembly_line = None
        self.architecture = None

        if name_restriction == None:
            self.pet_name = pn.generate_random_name()
        else:
            self.pet_name = pn.generate_random_name(name_restriction)


    def __lt__(self, other):
        """
        Override for less-than (<) comparison. Used for sorting individuals based on their 
        individual fitness. Whether the sorting happens in ascending or descending order is 
        determined by {self.optimization_type}.

        Example 1 - self.optimization_type=="max":
            self.fitness = 1.0
            other.fitness = 0.0
            --> check: 1.0 < 0.0 --> returns False
        Example 2 - self.optimization_type=="min":
            self.fitness = 1.0
            other.fitness = 0.0
            --> check: 0.0 < 1.0 --> returns True
        
        Args:
            other (_Individual): Another individual to compare to.

        Returns:
            boolean: 
                for {self.optimization_type=="min"}: True if other's fitness is lower than this one's, False otherwise.
                for {self.optimization_type=="max"}: True if this individual's fitness is lower than the other's, False otherwise.
        """
        if self.optimization_type=="max":
            return self.fitness < other.fitness
        else:
            return other.fitness < self.fitness


    def __str__(self, verbosity='vv'):
        """
        Returns a formatted string representation of the _Individual object.

        Args:
            verbosity (str, optional): The level of detail in the output. Defaults to 'vv'.

        Returns:
            str: A formatted string representation.
        """
        best_mol = None
        if len(self.molecules) > 0:
            best_mol = self.molecules[0].smiles_string

        ind = "Individual '%s':\n" % self.pet_name
        gen = "  Genotype:  %s\n" % self.genome
        gen_len = "  Genotype length:  %s\n" % len(self.genome)
        cod = "  Used codons:  %s\n" % self.used_codons
        phe = "  Phenotype: %s\n" % self.phenotype
        mol_count = "  Generated molecules: %s\n" % len(self.molecules)
        mol_best = "  Best molecule (SMILES): %s\n" % best_mol
        fit = "  Fitness:   %s\n" % self.fitness
        arch = "\nPKS Architecture: %s\n" % self.architecture

        if verbosity == 'v':
            # return the most compact values (useful when analyzing all inds and not just the best)
            return ind + gen_len + phe + mol_best + fit
        elif verbosity == 'vv':
            # default output: return everything except architecture (useful for best ind of each gen)           
            return ind + gen + gen_len + cod + phe  + mol_count + mol_best + fit
        elif verbosity == 'vvv':
            # return everything (useful for overall best individual in whole simulation)  
            return ind + gen + gen_len + cod + phe + mol_count + mol_best + fit + arch
        else:
            # return the default verbosity level (='vv') output
            return ind + gen + gen_len + cod + phe  + mol_count + mol_best + fit


    def sort_molecules(self):
        """
        This functions sorts the molecules of the Individual from best to worst based
        on the molecules fitness. What is considered best depends on the optimization_type
        of the individual.
        For maximization (optimization_type ='max') the list is ordered in descending order, otherwise
        the list will be sorted in ascending order.
        """
        # Sorting really only makes sense when there are at least two items
        if len(self.molecules) > 1:
            if self.optimization_type=="max":
                self.molecules.sort(key=lambda mol: mol.fitness, reverse=True) 
            else:
                self.molecules.sort(key=lambda mol: mol.fitness) 


    def build_pks_architecture(self):
        """
        Determine the modules and domains from the phenotype string.

        Returns:
            PKSArchitecture: The PKS architecture object representing the modules and domains.
        """

        if self.architecture:
            return self.architecture

        if not self.phenotype:
            # No phenotype? Then no architecture!
            return None
        
        module_list = []
        module_count = 0
        domain_count = 0

        # Replace dashes within parentheses with a temporary character, since inside those parantheses could be molecules that contain dashes
        # and we later want to only split at dashes, which are representing boundaries between domains
        input_string = self.phenotype
        temp_replaced_string = re.sub(r'AT\((.*?)\)', lambda match: match.group().replace("-", "~"), input_string)

        modules = temp_replaced_string.split("--")
        module_count = len(modules)
        
        for index, module in enumerate(modules):
            if index == 0:
                module_name = "Loading Module"
            else:
                module_name = f"{index}. Module"
            module_domains = module.split("-")

            # Replace the temporary character back to dashes within molecule names
            module_domains = [domain.replace("~", "-") for domain in module_domains]

            domain_count += len(module_domains)
            module_list.append([module_name, module_domains])

        self.architecture = PKSArchitecture(structure=module_list, module_count=module_count, domain_count=domain_count)  
        return self.architecture 