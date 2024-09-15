import json

from grammatical_evolution.model._individual import _Individual

class IndividualEncoder(json.JSONEncoder):
    """
    Class for encoding instances of class `grammatical_evolution.model._Individual` into JSON
    objects.
    """
    def default(self, obj: _Individual):
        if isinstance(obj, _Individual):
            return {
                'fitness': obj.fitness,
                'genome': obj.genome,
                'phenotype': obj.phenotype,
                'used_codons': obj.used_codons,
                'molecules': [
                    {
                        'smiles_string': molecule.smiles_string,
                        'fitness': molecule.fitness
                    }
                    for molecule in obj.molecules
                ]
                
            }
        return super().default(obj)