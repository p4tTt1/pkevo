# Note: Based on Robert Haas' panakeias_garden package (file: _enzyme_chemistry.py)
from external_tools.mod_importer import mod

class SubstrategyGenerator:
    """
    This class defines functions to generate and manage substrategies for chemical reactions (which are stored as .gml files).
    Substrategies are used for constructing and deconstructing (Note: Deconstruction will probably not be implemented in the 
    first release version of PKevo) Polyketide molecules. These substrategies are later put together into a single strategy 
    based on the used domains and substrates of the PKS evolved by the GrammaticalEvolution component.
    """
    def __init__(self, molecules, reactions):
        self.molecules = molecules
        self.reactions = reactions


    # Define a couple of shortcuts
    def _sequence(self, steps):
        return mod.DGStrat.makeSequence(steps)


    def _parallel(self, reactions):
        return mod.DGStrat.makeParallel(reactions)


    def _rule(self, reaction):
        return mod.DGStrat.makeRule(reaction)


    def _aromatization(self):
            return self._parallel([
                self._rule(self.reactions['fwd']['aromatization']['rs6_3db']),
                self._rule(self.reactions['fwd']['aromatization']['rs6_3ca']),
                self._rule(self.reactions['fwd']['aromatization']['rs6_2db_1ab']),
                self._rule(self.reactions['fwd']['aromatization']['rs6_1db_2ab']),
                self._rule(self.reactions['fwd']['aromatization']['rs6_2ca_1hy']),
                self._rule(self.reactions['fwd']['aromatization']['rs6_2ca_1db']),
                self._rule(self.reactions['fwd']['aromatization']['rs6_2ca_1ab']),
                self._rule(self.reactions['fwd']['aromatization']['rs6_1ca_2ab']),
                self._rule(self.reactions['fwd']['aromatization']['rs6_1ca_1db_1ab_ori1']),
                self._rule(self.reactions['fwd']['aromatization']['rs6_1ca_1db_1ab_ori2']),
                self._rule(self.reactions['fwd']['aromatization']['rs6_1ca_1db_1ab_ori3']),
                self._rule(self.reactions['fwd']['aromatization']['rs6_1ca_1hy_1ab_ori1']),
                self._rule(self.reactions['fwd']['aromatization']['rs6_1ca_1hy_1ab_ori2']),
                self._rule(self.reactions['fwd']['aromatization']['rs6_1ca_1hy_1db_ori1']),
                self._rule(self.reactions['fwd']['aromatization']['rs6_1ca_1hy_1db_ori2']),
                self._rule(self.reactions['fwd']['aromatization']['rs6_1ca_2db_ori1']),
                self._rule(self.reactions['fwd']['aromatization']['rs6_1ca_2db_ori2']),
                self._rule(self.reactions['fwd']['aromatization']['rs6_1hy_2ab']),
                self._rule(self.reactions['fwd']['aromatization']['rs6_1hy_2db_ori1']),
                self._rule(self.reactions['fwd']['aromatization']['rs6_1hy_2db_ori2']),
                self._rule(self.reactions['fwd']['aromatization']['rs6_1hy_1db_1ab_ori1']),
                self._rule(self.reactions['fwd']['aromatization']['rs6_1hy_1db_1ab_ori2']),
                self._rule(self.reactions['fwd']['aromatization']['rs6_1hy_1db_1ab_ori3']),
            ])


    def _beta_processing_preserving(self):
        return self._sequence([
            self._parallel([
                self._rule(self.reactions['fwd']['main_chain']['beta_carbonyl_reduction']),
                mod.filterSubset(True)
            ]),
            self._parallel([
                self._rule(self.reactions['fwd']['main_chain']['beta_hydroxyl_dehydration']),
                mod.filterSubset(True)
            ]),
            self._parallel([
                self._rule(self.reactions['fwd']['main_chain']['enoyl_reduction']),
                mod.filterSubset(True)
            ])
        ])


    def _beta_processing_nonpreserving(self):
        return (
            self.reactions['fwd']['main_chain']['beta_carbonyl_reduction']
            >> self.reactions['fwd']['main_chain']['beta_hydroxyl_dehydration']
            >> self.reactions['fwd']['main_chain']['enoyl_reduction']
        )


    def _cyclization(self):
        return self._parallel([
            self._rule(self.reactions['fwd']['cyclization']['aldol_addition_rs6_nc']),
            self._rule(self.reactions['fwd']['cyclization']['aldol_addition_rs6_dc']),
            self._rule(self.reactions['fwd']['cyclization']['aldol_condensation_rs6_nc']),
            self._rule(self.reactions['fwd']['cyclization']['aldol_condensation_rs6_dc']),
            self._rule(self.reactions['fwd']['cyclization']['decarboxylative_aldol_condensation_rs6']),
            #self._rule(self.reactions['fwd']['cyclization']['claisen_condensation_rs6_nc']),
            #self._rule(self.reactions['fwd']['cyclization']['claisen_condensation_rs6_dc']),
            #self._rule(self.reactions['fwd']['cyclization']['lactonization_aromatization_rs6']),
        ])

    
    def _post_pks(self):
        return self._parallel([
            self._rule(self.reactions['fwd']['post_pks']['epoxidation']),
            self._rule(self.reactions['fwd']['post_pks']['aliphatic_oxidation']),
            self._rule(self.reactions['fwd']['post_pks']['aromatic_oxidation']),
            self._rule(self.reactions['fwd']['post_pks']['o_glycosylation']),
            self._rule(self.reactions['fwd']['post_pks']['o_methylation']),
            self._rule(self.reactions['fwd']['post_pks']['o_carbamoylation']),
        ])

    
    def _release_hydroxyl_group_water(self):
        # (Robert) TODO: faster? only monomorphism if same atom count (or part of mod function)?
        return mod.leftPredicate[lambda derivation:
            any(g.isomorphism(self.molecules['basic_units']['Water']) for g in derivation.left)](
                self.reactions['fwd']['release']['hydroxyl'])


    def _release_hydroxyl_group_intramolecular(self):
        return mod.leftPredicate[lambda derivation: len(derivation.left) == 1](
            self.reactions['fwd']['release']['hydroxyl'])

    
    def _release_amino_group_intramolecular(self):
        return mod.leftPredicate[lambda aDerivation: len(aDerivation.left) == 1](
            self.reactions['fwd']['release']['amino'])

    
    def _release_claisen(self):
        return self._parallel([
            self._rule(self.reactions['fwd']['release']['claisen_condensation_rs6_nc']),
            self._rule(self.reactions['fwd']['release']['claisen_condensation_rs6_dc']),
        ])


    def generate_substrategies(self):
        """
        Generate the defined Substrategies.

        Returns:
            dict: A dictionary containing forward and reverse substrategies.
        """
        substrategies = {
            'fwd': dict(), # forward strategies are for constructing the molecules (Domain String -> Molecule)
            'rev': dict(), # reverse strategies are for deconstructing the molecules (Molecule -> Domain String)
        }

        # TODO: Still a bit messy. Some "substrategies" are generated directly inside the PK Strategy generator
        substrategies['fwd']['aromatization'] = self._aromatization()
        # substrategies['fwd']['beta_processing_nonpreserving'] = self._beta_processing_nonpreserving()
        # substrategies['fwd']['beta_processing_preserving'] = self._beta_processing_preserving()
        substrategies['fwd']['cyclization'] = self._cyclization()
        substrategies['fwd']['post_pks'] = self._post_pks()
        substrategies['fwd']['release_amino_group_intramolecular'] = self._release_amino_group_intramolecular()
        substrategies['fwd']['release_claisen'] = self._release_claisen()
        substrategies['fwd']['release_hydroxyl_group_intramolecular'] = self._release_hydroxyl_group_intramolecular()
        substrategies['fwd']['release_hydroxyl_group_water'] = self._release_hydroxyl_group_water()

        return substrategies