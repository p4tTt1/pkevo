from dataclasses import dataclass

@dataclass(frozen=False)
class Molecule:
    """
    A simple lass representing chemical molecules by a SMILES string instead of
    a more complex graph, like MØD uses to do. In addition, the score received by the
    GE evalution is stored for the molecule as well. This score can potentially change,
    the SMILES representation should not (hence the leading underscore to indicate this. 
    Python unfortunately doesn't offer a fool-proof way of forbiding individual attributes 
    to be READ-ONLY). 
    """
    _smiles_string: str
    fitness: float

    def __init__(self, smiles_string: str, fitness: float):
        self._smiles_string = smiles_string
        self.fitness = fitness

    @property
    def smiles_string(self):
        return self._smiles_string

    def __str__(self):
        return f"{self.smiles_string}, Fitness score: {self.fitness}"