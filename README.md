# parsePrIMe

Run `parse_prime.py` with the appropriate arguments to parse all reactions in a PrIMe kinetics database and convert them to RMG reactions. Structural information is only available in InChI format or through CAS numbers in the database; for InChI, direct conversion is attempted using RMG, and for CAS, CIRpy is used to pull the structural information. For the current version of the database, this means that approximately 600 species cannot be parsed.
