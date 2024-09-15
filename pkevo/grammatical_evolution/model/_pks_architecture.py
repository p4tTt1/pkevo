from dataclasses import dataclass

@dataclass(frozen=True)
class PKSArchitecture:
    """
    Represents the architecture of Polyketide Synthases (PKS).

    Attributes:
        structure (list): A hierarchical structure describing the arrangement of modules and domains in the PKS.
        module_count (int): The total number of modules in the PKS.
        domain_count (int): The total number of domains in the PKS.

    Methods:
        __str__():
            Returns a formatted string representation of the PKSArchitecture object.

        format_structure(substructure, level):
            Formats the hierarchical structure of the PKS for display.

    Note:
        The PKSArchitecture class is intended to provide a structured representation of PKS architecture.
        It is typically used for visualization and analysis of PKS structures.
    """
    structure: list()
    module_count: int
    domain_count: int


    def __str__(self):
        """
        Returns a formatted string representation of the PKSArchitecture object.
        """
        struct = ["\n.",
                  *self.format_structure(self.structure, 0)]
        mods = f"\n\nNumber of Modules: {self.module_count}\n"
        doms = f"Number of Domains: {self.domain_count}\n"
        return '\n'.join(struct) + mods + doms


    def format_structure(self, substructure, level):
        """
        Formats the hierarchical structure of the PKS for display.

        Args:
            substructure (list): The substructure of the PKS.
            level (int): The current level of indentation in the hierarchical structure.

        Returns:
            list: A list of strings representing the formatted hierarchical structure.
        """
        result = []
        for index, item in enumerate(substructure):
            if isinstance(item, list):
                connector = "└── " if index == len(substructure) - 1 else "├── "
                result.append("│ " * level + connector + item[0])
                result.extend(self.format_structure(item[1], level + 1))
            else:
                connector = "└── " if index == len(substructure) - 1 else "├── "
                result.append("│ " * level + connector + item)
        return result