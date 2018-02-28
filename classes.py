#!/usr/bin/env python
# -*- coding:utf-8 -*-

import cirpy

from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.reaction import Reaction


class ConversionError(Exception):
    """
    An exception class for errors during name to structure conversion.
    """
    pass


class StoichiometryError(Exception):
    """
    An exception class for invalid reaction stoichiometries.
    """


class PrIMeSpecies(object):
    def __init__(self, prime_id, inchi=None, cas=None, names=None, composition=None):
        """
        composition should be a dictionary of atomic symbols and counts.
        """
        self.prime_id = prime_id
        self.inchi = inchi
        self.cas = cas
        self.names = names
        self.composition = composition

    def __str__(self):
        return self.prime_id

    def __repr__(self):
        return self.__class__.__name__ + '(' + str(self) + ')'

    def get_rmg_species(self):
        """
        Convert identifiers from PrIMe database to RMG molecule with
        structural information. Try InChI first, then CAS look-up,
        and finally try to parse species names.
        """
        if self.inchi is not None:
            return self.parse_inchi()

        cas_error = False
        if self.cas is not None:
            try:
                return self.parse_cas()
            except ConversionError:
                if self.names is not None:
                    cas_error = True
                else:
                    raise

        if self.names is not None:
            try:
                return self.parse_names()
            except ConversionError:
                if cas_error:
                    raise ConversionError('Could not resolve CAS and name for species {}.'.format(self.prime_id))
                else:
                    raise

        raise ConversionError(
            'InChI, CAS, or name required for structure conversion of species {}.'.format(self.prime_id)
        )

    def parse_inchi(self):
        rmg_species = Species()
        rmg_species.molecule = [Molecule().fromInChI(self.inchi, backend='rdkit')]
        return rmg_species
            
    def parse_cas(self):
        smiles = cirpy.resolve(self.cas, 'smiles', ['cas_number'])
        if smiles is None:
            raise ConversionError('Could not resolve CAS number for species {}.'.format(self.prime_id))
        else:
            return Species().fromSMILES(smiles)

    def parse_names(self):
        for name in self.names:
            smiles = cirpy.resolve(name, 'smiles', ['name_by_opsin'])
            if smiles is not None:
                return Species().fromSMILES(smiles)
        else:
            for name in self.names:
                smiles = cirpy.resolve(name, 'smiles', ['name_by_cir'])
                if smiles is not None:
                    return Species().fromSMILES(smiles)
            else:
                raise ConversionError('Could not resolve name for species {}.'.format(self.prime_id))


class PrIMeReaction(object):
    def __init__(self, prime_id, species_coeffs, species_dict):
        """
        species_coeffs should be a list of tuples each containing
        the prime_id of a species and its stoichiometric coefficient.
        The dictionary of PrIMeSpecies is required to match the
        structure information to the reaction.
        """
        self.prime_id = prime_id
        self.reactants = None
        self.products = None
        self.rmg_reactants = None
        self.rmg_products = None
        self.kinetics = None
        self._direction_matched = False

        self._parse_species_coeffs(species_coeffs, species_dict)
        self._check_stoichiometry()

    def __str__(self):
        return self.prime_id

    def __repr__(self):
        return self.__class__.__name__ + '(' + str(self) + ')'

    def _parse_species_coeffs(self, species_coeffs, species_dict):
        """
        Create reactant and product lists using the species PrIMe IDs
        and dictionary of PrIMeSpecies.
        """
        self.reactants = []
        self.products = []
        
        for species_id, coeff in species_coeffs:
            try:
                if coeff < 0:
                    self.reactants.extend(abs(coeff)*[species_dict[species_id]])
                elif coeff > 0:
                    self.products.extend(coeff*[species_dict[species_id]])
            except KeyError:
                raise KeyError('Species {} not found in species dictionary.'.format(species_id))

    def _check_stoichiometry(self):
        reactant_counts = {}
        product_counts = {}

        for reactant in self.reactants:
            if reactant.composition is None:
                raise StoichiometryError(
                    'Species {} is missing composition in reaction {}.'.format(reactant.prime_id, self.prime_id)
                )
            for atom, count in reactant.composition.iteritems():
                if atom in reactant_counts:
                    reactant_counts[atom] += count
                else:
                    reactant_counts[atom] = count
        for product in self.products:
            if product.composition is None:
                raise StoichiometryError(
                    'Species {} is missing composition in reaction {}.'.format(product.prime_id, self.prime_id)
                )
            for atom, count in product.composition.iteritems():
                if atom in product_counts:
                    product_counts[atom] += count
                else:
                    product_counts[atom] = count

        if not reactant_counts == product_counts:
            raise StoichiometryError('Reaction {} has invalid stoichiometry.'.format(self.prime_id))

    def match_direction(self):
        """
        Match the direction of reactants and products to match the
        kinetics expression.
        """
        if self.reactants is None or self.products is None:
            raise Exception('Missing reactants or products for reaction {}.'.format(self.prime_id))

        if self.kinetics is not None:
            if self.kinetics.direction == 'forward':
                self._direction_matched = True
            elif self.kinetics.direction == 'reverse':
                new_products = list(self.reactants)
                self.reactants = list(self.products)
                self.products = new_products

                if self.rmg_reactants is not None and self.rmg_products is not None:
                    new_rmg_products = list(self.rmg_reactants)
                    self.rmg_reactants = list(self.rmg_products)
                    self.rmg_products = new_rmg_products

                self._direction_matched = True
            else:
                raise Exception(
                    'Invalid kinetics direction "{}" in reaction {}.'.format(self.kinetics.direction, self.prime_id)
                )
        else:
            raise Exception('No kinetics set for reaction {}.'.format(self.prime_id))

    def get_rmg_species_from_dict(self, rmg_species_dict):
        """
        Given a dictionary of RMG species with the same PrIMe ID
        keys, extract the species corresponding to this reaction.
        """
        self.rmg_reactants = []
        self.rmg_products = []

        for spc in self.reactants:
            self.rmg_reactants.append(rmg_species_dict[spc.prime_id])
        for spc in self.products:
            self.rmg_products.append(rmg_species_dict[spc.prime_id])

    def get_rmg_reaction(self):
        """
        Convert PrIMeReaction to an RMG reaction.
        """
        if not all((self.reactants, self.products, self.kinetics)):
            raise ConversionError('Missing reactants, products, or kinetics for reaction {}.'.format(self.prime_id))
        if not self._direction_matched:
            self.match_direction()
        
        reaction = Reaction()
        if self.rmg_reactants is None or self.rmg_products is None:
            reaction.reactants = [s.get_rmg_species() for s in self.reactants]
            reaction.products = [s.get_rmg_species() for s in self.products]
        else:
            reaction.reactants = self.rmg_reactants
            reaction.products = self.rmg_products
        reaction.kinetics = self.kinetics.expression

        return reaction


class PrIMeKinetics(object):
    def __init__(self, prime_id,
                 expression=None, has_n=False, rate_law_type=None, direction=None, year=None, links=None):
        self.prime_id = prime_id
        self.expression = expression
        self.has_n = has_n
        self.rate_law_type = rate_law_type
        self.direction = direction
        self.year = year
        self.links = links
