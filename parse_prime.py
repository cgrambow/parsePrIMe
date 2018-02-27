#!/usr/bin/env python
# -*- coding:utf-8 -*-

import argparse
import glob
import os
import re
import warnings
import xml.etree.cElementTree as ET

from rmgpy.kinetics import Arrhenius

from classes import PrIMeSpecies, PrIMeReaction, PrIMeKinetics, StoichiometryError
from util import custom_warning

warnings.formatwarning = custom_warning

# Compile regular expression searches
nsprog = re.compile(r'\{.*\}')          # Matches namespaces
yearprog = re.compile(r'(19|20)\d{2}')  # Matches years

# Other settings/constants
valid_params = {'a', 'e', 'n'}


class KineticsError(Exception):
    """
    An exception class for errors during kinetics parsing.
    """
    pass


def main():
    parser = argparse.ArgumentParser(
        description='Parse the PrIMe database.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('depository', type=str, nargs='?', default='depository', metavar='PATH',
                        help='Path to PrIMe depository')
    args = parser.parse_args()
    depository = args.depository

    species_dict = parse_species(depository)
    reactions_dict = parse_reactions(depository, species_dict)
    nrxn = len(reactions_dict)
    reactions_dict = get_kinetics(depository, reactions_dict)

    print('Number of valid reactions: {}'.format(nrxn))
    print('Number of valid reactions with kinetics: {}'.format(len(reactions_dict)))


def parse_species(depository):
    """
    Loop through all species XMLs and extract chemical identifiers
    and chemical compositions. Save the information as attributes of
    PrIMeSpecies and group in a dictionary with PrIMe ID keys.
    """
    species_dir = os.path.join(depository, 'species', 'catalog')
    species_dict = {}

    for xml in glob.iglob(os.path.join(species_dir, 's*.xml')):
        try:
            spc = parse_single_species(xml)
        except ET.ParseError:
            warnings.warn('Skipped {}: invalid XML.'.format(xml))
            continue
        species_dict[spc.prime_id] = spc

    return species_dict


def parse_reactions(depository, species_dict):
    """
    Loop through all reaction XMLs and extract species participating
    in the reaction. Save the information in PrIMeReaction objects
    and group in a dictionary with PrIMe ID keys.
    """
    reactions_dir = os.path.join(depository, 'reactions', 'catalog')
    reactions_dict = {}

    for xml in glob.iglob(os.path.join(reactions_dir, 'r*.xml')):
        try:
            rxn = parse_single_reaction(xml, species_dict)
        except ET.ParseError:
            warnings.warn('Skipped {}: invalid XML.'.format(xml))
            continue
        except StoichiometryError as e:
            warnings.warn('Skipped {}: {}'.format(xml, e))
        reactions_dict[rxn.prime_id] = rxn

    return reactions_dict


def get_kinetics(depository, reactions_dict):
    """
    Parse XML files containing kinetics data for each reaction.
    Only returns reactions with at least Arrhenius kinetics.
    """
    data_parent_dir = os.path.join(depository, 'reactions', 'data')
    reactions_dict_out = {}

    for prime_id, rxn in reactions_dict.iteritems():
        data_dir = os.path.join(data_parent_dir, prime_id)
        kinetics_list = []

        for xml in glob.iglob(os.path.join(data_dir, 'rk*.xml')):
            try:
                prime_kinetics = parse_kinetics(xml, rxn)
                # Use newest reference, preferably in forward direction, also first try to use ones with n
                # If rate_law_type is sum, ignore links
            except ET.ParseError:
                warnings.warn('Skipped {}: invalid XML.'.format(xml))
                continue
            except KineticsError as e:
                warnings.warn('Skipped {}: {}'.format(xml, e))
            else:
                kinetics_list.append(prime_kinetics)

        # Determine which kinetics expression to keep
        if kinetics_list:
            # Remove kinetics corresponding to a "sum" rate law type
            sum_types = {k for k in kinetics_list if k.rate_law_type == 'sum'}
            links = {kl for kl in kinetics_list for k in sum_types if kl.prime_id in kinetics_list[k].links}
            kinetics_list_trimmed = [k for k in kinetics_list if k not in sum_types and k not in links]

            if kinetics_list_trimmed:
                kinetics_list_trimmed.sort(key=lambda kin: kin.year, reverse=True)  # Use newest data first
                kinetics_with_n = [k for k in kinetics_list_trimmed if k.has_n]  # Prefer to use data that has n

                if kinetics_with_n:
                    kinetics_list_final = kinetics_with_n
                else:
                    kinetics_list_final = kinetics_list_trimmed

                # Prefer kinetics that have a year and are in the forward direction
                for k in kinetics_list_final:
                    if k.year is not None and k.direction == 'forward':
                        chosen_kinetics = k
                        break
                else:
                    for k in kinetics_list_final:
                        if k.year is not None:
                            chosen_kinetics = k
                            break
                    else:
                        chosen_kinetics = kinetics_list_final[0]

                rxn.kinetics = chosen_kinetics
                reactions_dict_out[prime_id] = rxn

    return reactions_dict_out


def parse_single_species(xml):
    """
    Parse the information from a single species XML.
    """
    tree = ET.parse(xml)
    root = tree.getroot()
    prime_id = root.attrib['primeID']  # Extract ID

    # Extract namespace
    ns = nsprog.match(root.tag).group(0)

    # Find CAS number and InChI if available
    identifiers = root.find(ns + 'chemicalIdentifier')
    cas = None
    inchi = None
    names = []
    for identifier in identifiers:
        try:
            identifier_type = identifier.attrib['type'].lower()
        except KeyError:
            names.append(identifier.text)
            continue
        if identifier_type == 'casregistrynumber':
            cas = identifier.text
        elif identifier_type == 'inchi':
            inchi = identifier.text
    if not names:
        names = None

    # Extract chemical composition
    atoms = root.find(ns + 'chemicalComposition')
    composition = {}
    for atom in atoms:
        if atom.tag == ns + 'atom':
            composition[atom.attrib['symbol'].upper()] = int(atom.text)
    if not composition:
        composition = None

    return PrIMeSpecies(prime_id, inchi=inchi, cas=cas, names=names, composition=composition)


def parse_single_reaction(xml, species_dict):
    """
    Parse the information from a single reaction XML.
    """
    tree = ET.parse(xml)
    root = tree.getroot()
    prime_id = root.attrib['primeID']  # Extract ID

    # Extract namespace
    ns = nsprog.match(root.tag).group(0)

    # Extract species
    species = root.find(ns + 'reactants')
    species_coeffs = []
    for spc in species:
        species_coeffs.append((spc.attrib['primeID'], int(spc.text)))

    return PrIMeReaction(prime_id, species_coeffs, species_dict)


def parse_kinetics(xml, rxn, pdep=False):
    """
    Parse the information from a single kinetics XML. If a P-dep
    expression is encountered and pdep is set to False, a
    KineticsError will be raised.
    """
    tree = ET.parse(xml)
    root = tree.getroot()
    prime_id = root.attrib['primeID']  # Extract ID

    # Extract namespace
    ns = nsprog.match(root.tag).group(0)

    # Make sure kinetics exist
    rate_coeff = root.find(ns + 'rateCoefficient')
    if rate_coeff is None:
        raise KineticsError('No rate coefficient found.')

    # Extract year
    ref = root.find(ns + 'bibliographyLink')
    year_match = yearprog.search(ref.attrib['preferredKey'])
    if year_match is None:
        year = None
        warnings.warn('No year in {}.'.format(xml))
    else:
        year = int(year_match.group(0))

    # Extract rate law type
    rate_law_type = root.attrib['rateLawType']  # Check for key error?
    if rate_law_type == 'mass action':
        pass
    elif rate_law_type == 'sum':
        # Save the links to the kinetics that are meant to be summed up
        reaction_link = root.find(ns + 'reactionLink')
        links = {rate_link.attrib['primeID'] for rate_link in reaction_link}
        return PrIMeKinetics(prime_id, rate_law_type='sum', links=links)
    elif rate_law_type in ('third body', 'unimolecular', 'chemical activation'):
        if pdep:
            raise NotImplementedError('P-dep has not been implemented for "{}" type.'.format(rate_law_type))
        else:
            raise KineticsError('P-dep kinetics: {}.'.format(rate_law_type))
    else:
        raise Exception('Unknown rate law type: {}'.format(rate_law_type))

    # Extract kinetics
    direction = rate_coeff.attrib['direction'].lower()

    expression = rate_coeff.find(ns + 'expression')
    form = expression.attrib['form'].lower()
    if form == 'constant':
        raise KineticsError('Missing activation energy.')
    elif form != 'arrhenius':
        raise NotImplementedError('Kinetics form "{}" has not been implemented yet.'.format(form))

    a = e = None
    n = 0.0
    has_n = False
    for param in expression:
        if param.tag == ns + 'parameter':
            if param.attrib['name'].lower() not in valid_params:
                raise Exception('Unknown parameter: {}'.format(param.attrib['name']))
            if param.attrib['name'].lower() == 'a':
                units = _format_units(param.attrib['units'], rxn, direction)
                a = (float(param.find(ns + 'value').text), units)
            elif param.attrib['name'].lower() == 'e':
                units = _format_units(param.attrib['units'], rxn, direction)
                e = (float(param.find(ns + 'value').text), units)
            elif param.attrib['name'].lower() == 'n':
                n = float(param.find(ns + 'value').text)
                has_n = True

    if a is None:
        raise KineticsError('Missing prefactor.')
    if e is None:
        raise KineticsError('Missing activation energy.')

    # If A has units involving mol and is very small, then it's probably molecules
    if 'mol' in a[1] and a[0] < 1.0:
        a[1].replace('mol', 'molecule')
        warnings.warn('Replaced "mol" by "molecule" in units for {}.'.format(xml))

    return PrIMeKinetics(prime_id,
                         expression=Arrhenius(A=a, n=n, Ea=e),
                         has_n=has_n,
                         rate_law_type=rate_law_type,
                         direction=direction,
                         year=year)


def _format_units(units, rxn, direction):
    """
    Given a string like "cm3,mol,s,K" for the units of a rate
    coefficient, format it to have the correct form using the
    reaction molecularity. Can also pass activation energy units to
    check if they are missing.
    """
    if units == 's':
        return 's^-1'
    elif units == 'cm3,mol,s,K' or units == 'cm3,mol,s':
        if direction == 'forward':
            if len(rxn.reactants) == 2:
                return 'cm^3/(mol*s)'
            else:
                # NOTE: Changed this to kinetics error to avoid possibly incorrect kinetics for now
                raise KineticsError(
                    '{} reactants have not been implemented yet for units {}.'.format(len(rxn.reactants), units)
                )
        elif direction == 'reverse':
            if len(rxn.products) == 2:
                return 'cm^3/(mol*s)'
            else:
                # NOTE: Changed this to kinetics error to avoid possibly incorrect kinetics for now
                raise KineticsError(
                    '{} products have not been implemented yet for units {}.'.format(len(rxn.products), units)
                )
    elif units == 'cm3,molecule,s,K' or units == 'cm3,molecule,s':
        if direction == 'forward':
            if len(rxn.reactants) == 2:
                return 'cm^3/(molecule*s)'
            else:
                # NOTE: Changed this to kinetics error to avoid possibly incorrect kinetics for now
                raise KineticsError(
                    '{} reactants have not been implemented yet for units {}.'.format(len(rxn.reactants), units)
                )
        elif direction == 'reverse':
            if len(rxn.products) == 2:
                return 'cm^3/(molecule*s)'
            else:
                # NOTE: Changed this to kinetics error to avoid possibly incorrect kinetics for now
                raise KineticsError(
                    '{} products have not been implemented yet for units {}.'.format(len(rxn.products), units)
                )
    elif units == 'cm6,mol,s':
        if direction == 'forward':
            if len(rxn.reactants) == 3:
                return 'cm^6/(mol^2*s)'
            else:
                raise NotImplementedError(
                    '{} reactants have not been implemented yet for units {}.'.format(len(rxn.reactants), units)
                )
        elif direction == 'reverse':
            if len(rxn.products) == 3:
                return 'cm^6/(mol^2*s)'
            else:
                raise NotImplementedError(
                    '{} products have not been implemented yet for units {}.'.format(len(rxn.products), units)
                )
    elif units == 'K' or units == 'cal/mol' or units == 'kcal/mol':
        return units
    elif units == '':
        raise KineticsError('Missing units.')
    else:
        raise NotImplementedError('Units {} have not been implemented yet.'.format(units))


if __name__ == '__main__':
    main()
