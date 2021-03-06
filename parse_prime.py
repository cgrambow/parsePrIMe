#!/usr/bin/env python
# -*- coding:utf-8 -*-

from __future__ import print_function

import argparse
import glob
import os
import re
import sys
import urllib2
import warnings
import xml.etree.cElementTree as ET

from rmgpy.kinetics import Arrhenius
from rmgpy.exceptions import AtomTypeError

from classes import PrIMeSpecies, PrIMeReaction, PrIMeKinetics, ConversionError, StoichiometryError
from util import custom_warning, pickle_dump, pickle_load

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
    # Set up command line arguments
    parser = argparse.ArgumentParser(
        description='Parse the PrIMe database.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', '--outdir', type=str, metavar='DIR',
                        help='Directory to save pickled dictionaries and lists in')
    parser.add_argument('depository', type=str, nargs='?', default='depository', metavar='PATH',
                        help='Path to PrIMe depository')
    args = parser.parse_args()
    out_dir = args.outdir
    depository = args.depository

    prime_species_dict = prime_reactions_dict = prime_species_in_reactions_dict = rmg_species_dict = None
    prime_species_path = os.path.join(out_dir, 'prime_species_dict.pickle')
    prime_reactions_path = os.path.join(out_dir, 'prime_reactions_dict.pickle')
    prime_species_in_reactions_path = os.path.join(out_dir, 'prime_species_in_reactions_dict.pickle')
    rmg_species_path = os.path.join(out_dir, 'rmg_species_dict.pickle')
    if out_dir is not None:
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        else:
            # Check if we can load existing dictionaries, which could save a lot of time
            print('Trying to load dictionaries...')
            try:
                prime_species_dict = pickle_load(prime_species_path)
            except IOError:
                print('Could not find prime_species_dict')
            else:
                print('Successfully loaded prime_species_dict')
            try:
                prime_reactions_dict = pickle_load(prime_reactions_path)
            except IOError:
                print('Could not find prime_reactions_dict')
            else:
                print('Successfully loaded prime_reactions_dict')
            try:
                prime_species_in_reactions_dict = pickle_load(prime_species_in_reactions_path)
            except IOError:
                print('Could not find prime_species_in_reactions_dict')
            else:
                print('Successfully loaded prime_species_in_reactions_dict')
            try:
                rmg_species_dict = pickle_load(rmg_species_path)
            except IOError:
                print('Could not find rmg_species_dict')
            else:
                print('Successfully loaded rmg_species_dict')

    if prime_species_dict is None:
        print('Parsing species...')
        prime_species_dict = parse_species(depository)

    parsed_reactions = False
    nrxn_pre_kinetics = None
    if prime_reactions_dict is None:
        print('Parsing reactions...')
        prime_reactions_dict = parse_reactions(depository, prime_species_dict)
        nrxn_pre_kinetics = len(prime_reactions_dict)
        print('Parsing kinetics...')
        prime_reactions_dict = get_kinetics(depository, prime_reactions_dict)
        parsed_reactions = True
        # Note: The reactions are not necessarily in the correct direction at this point.
        #       Have to run match_direction first.

    print('Number of valid PrIMe species: {}'.format(len(prime_species_dict)))
    if parsed_reactions:
        print('Number of valid PrIMe reactions: {}'.format(nrxn_pre_kinetics))
    print('Number of valid PrIMe reactions with kinetics: {}'.format(len(prime_reactions_dict)))

    if out_dir is not None:
        print('Saving PrIMe species and reactions dictionaries to {}'.format(out_dir))
        pickle_dump(prime_species_path, prime_species_dict)
        pickle_dump(prime_reactions_path, prime_reactions_dict)

    if prime_species_in_reactions_dict is None:
        print('Extracting species in reactions...')
        # Only convert species actually involved in reactions
        prime_species_in_reactions_dict = {}
        for rxn in prime_reactions_dict.itervalues():
            for spc in rxn.reactants:
                prime_species_in_reactions_dict[spc.prime_id] = spc
            for spc in rxn.products:
                prime_species_in_reactions_dict[spc.prime_id] = spc

    if out_dir is not None:
        print('Saving PrIMe species in reactions dictionary to {}'.format(out_dir))
        pickle_dump(prime_species_in_reactions_path, prime_species_in_reactions_dict)

    print('Converting species to RMG types...')
    if rmg_species_dict is None:
        rmg_species_dict = {}
    count_resolve_errors = 0
    for prime_id, spc in prime_species_in_reactions_dict.iteritems():
        # Don't bother converting if we already did so in a previous run
        if prime_id in rmg_species_dict:
            continue

        try:
            rmg_species_dict[prime_id] = spc.get_rmg_species()
        except ConversionError as e:
            count_resolve_errors += 1
            warnings.warn('Skipped {}: {}'.format(prime_id, e))
            continue
        except (ValueError, AttributeError, AtomTypeError) as e:
            warnings.warn('Skipped {}: {}'.format(prime_id, e))
            continue
        except KeyError as e:
            warnings.warn('Skipped {}: Atom type {} is not supported.'.format(prime_id, e))
            continue
        except urllib2.URLError as e:
            warnings.warn('URLError encountered for {}: {}, retrying...'.format(prime_id, e))
            try:
                rmg_species_dict[prime_id] = spc.get_rmg_species()
            except urllib2.URLError as e:
                warnings.warn('URLError encountered for {}: {}, retrying...'.format(prime_id, e))
                try:
                    rmg_species_dict[prime_id] = spc.get_rmg_species()
                except urllib2.URLError as e:
                    warnings.warn('Skipped {}: {}'.format(prime_id, e))
                    continue
        except Exception as e:
            if "Couldn't parse" in str(e):
                warnings.warn('Skipped {}: {}'.format(prime_id, e))
                continue
            else:
                print('Error encountered during conversion of species {}.'.format(prime_id), file=sys.stderr)
                # Save output regardless, so we don't have to do all the work again next time
                if out_dir is not None:
                    print('Saving RMG species dictionary to {}'.format(out_dir))
                    pickle_dump(rmg_species_path, rmg_species_dict)
                raise
        else:
            print('Converted {}.'.format(prime_id))

    print('Number of PrIMe species in reactions: {}'.format(len(prime_species_in_reactions_dict)))
    print('Number of RMG species in reactions: {}'.format(len(rmg_species_dict)))
    print('Number of CIRpy resolve errors: {}'.format(count_resolve_errors))

    if out_dir is not None:
        print('Saving RMG species dictionary to {}'.format(out_dir))
        pickle_dump(rmg_species_path, rmg_species_dict)

    print('Converting reactions to RMG types...')
    reactions = []
    for rxn in prime_reactions_dict.itervalues():
        try:
            rxn.get_rmg_species_from_dict(rmg_species_dict)
        except KeyError:
            continue
        else:
            reactions.append(rxn.get_rmg_reaction())

    print('Number of RMG reactions: {}'.format(len(reactions)))

    if out_dir is not None:
        print('Saving RMG reactions list to {}'.format(out_dir))
        pickle_dump(os.path.join(out_dir, 'reactions.pickle'), reactions)


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
        kinetics_dict = {}

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
                kinetics_dict[prime_kinetics.prime_id] = prime_kinetics

        # Determine which kinetics expression to keep
        if kinetics_dict:
            # Remove kinetics corresponding to a "sum" rate law type
            sum_types = {k for k in kinetics_dict.itervalues() if k.rate_law_type == 'sum'}
            links = {kl for kl in kinetics_dict.itervalues() for k in sum_types
                     if kl.prime_id in kinetics_dict[k.prime_id].links}
            kinetics_list_trimmed = [k for k in kinetics_dict.itervalues() if k not in sum_types and k not in links]

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

    # Extract rate law type
    try:
        rate_law_type = root.attrib['rateLawType']
    except KeyError:
        # Check later if there is a rate coefficient
        rate_law_type = None
    else:
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

    # Make sure kinetics exist
    rate_coeff = root.find(ns + 'rateCoefficient')
    if rate_coeff is None:
        raise KineticsError('No rate coefficient found.')

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
    # Smallness is assessed based on an "effective" A at 1000 K (i.e. A*T^n)
    if a[1] == 'cm^3/(mol*s)' and a[0] * 1000.0**n < 1.0:
        a[1].replace('mol', 'molecule')
        warnings.warn('Replaced "mol" by "molecule" in units for {}.'.format(xml))

    # Extract year
    ref = root.find(ns + 'bibliographyLink')
    year_match = yearprog.search(ref.attrib['preferredKey'])
    if year_match is None:
        year = None
        warnings.warn('No year in {}.'.format(xml))
    else:
        year = int(year_match.group(0))

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
