#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Little script to prepare data for kohonen learngin method.

Usage:
    prepare_kohonen.py pdb [-v --verbose] <ids_file>...
    prepare_kohonen.py dssp [-v --verbose] <pdb_file>...
    prepare_kohonen.py convert [-v --verbose -r --rm] <dssp_file>...
    prepare_kohonen.py angle
    prepare_kohonen.py (-h | --help | --version)

Arguments:
    <ids_file>  File with pdb's id one per line.
    <pdb_file>  File from Protein Data Bank.

Options:
    --version           Give the script version number.
    -h --help           Display the help.
    -v --verbose        Give more informations about script processus.
    -r --rm             Remove mutual.csv for convert command

"""

from __future__ import print_function

import os
import csv
import sys
import docopt
import requests
import subprocess

from time import strftime

# URL to download pdb file by their id
PDB_URL = "http://www.rcsb.org/pdb/files/{}.pdb"

# Directory to store downloaded pdbs
PDB_DIR = os.path.join(os.path.abspath("."), 'pdbs/')

# Location of the dssp executable
DSSP = "./dssp-2.0.4-linux-i386"

# Directory to store dssp file
DSSP_DIR = os.path.join(os.path.abspath("."), 'dssp/')

# Launcher of dssp_to_csv bash script
DSSP_TO_CSV = "./scripts/dssp_to_csv.sh"

# Location of mutual.csv
MUTUAL = "results/mutual.csv"

# Location of angles.csv
ANGLES = "results/angles.csv"

# Location of vectors.csv
VECTORS = "results/vectors.csv"


def is_file(files):
    """Check if the user give files as argument"""
    for ids_file in files:
        if not os.path.isfile(ids_file):
            print("{} isn't a file.".format(ids_file), file=sys.stderr)
            sys.exit()
    return True


def verify_file(files, ext='.txt'):
    """Check if the user give pdb files as argument"""
    for pdb_file in files:
        file_extension = os.path.splitext(pdb_file)[1]
        if not file_extension == ext or not is_file([pdb_file]):
            print("{} isn't a {} file.".format(pdb_file, ext[1:]),
                  file=sys.stderr)
            sys.exit()
    return True


def log_error(error_message):
    """Writes the occurred errors in a log file"""
    with open('log_error.txt', 'a') as f_log:
        f_log.write("[{}] ".format(strftime("%Y-%m-%d %H:%M:%S")))
        f_log.write(error_message + "\n")


def download_pdbs(args):
    """
    Download pdbs by their ids contained in files given by the user.
    PDB files are organized in directory named in function of the file of id
    used.
    """
    nb_errors = 0

    # Check if the pdbs' directory exists and create it if not
    if not os.path.isdir(PDB_DIR):
        os.mkdir(PDB_DIR)

    for ids_file in args['<ids_file>']:
        # Create separate directories for each file given to the script
        pdb_file_dir = os.path.join(PDB_DIR, os.path.splitext(ids_file)[0]+"/")
        if not os.path.isdir(pdb_file_dir):
            os.mkdir(pdb_file_dir)

        print("Getting ids from {} begins".format(ids_file))

        # Read each pdb's id
        with open(ids_file, 'r') as f_ids:
            for pdb_id in f_ids:
                # Case with blank line
                if pdb_id is None:
                    continue

                # Deletes the carriage return from the id
                pdb_id = pdb_id.rstrip()
                # Name of the downloaded pdb file
                file_name = pdb_file_dir + pdb_id + ".pdb"

                # Requests for the pdb file by id
                response = requests.get(PDB_URL.format(pdb_id))

                # unprocessed request
                if not response.ok:
                    nb_errors += 1
                    message = "Error: Could not download the pdb " +\
                              "file: {}.pdb".format(pdb_id)
                    print(message, file=sys.stderr)
                    log_error(message)
                    continue

                # Wrong id
                if "requested file is not available." in response.text:
                    nb_errors += 1
                    message = "The id: {} doesn't seem to exist.".format(pdb_id)
                    print(message, file=sys.stderr)
                    log_error(message)
                    continue

                # Writes the pdb file
                with open(file_name, 'w') as f_pdbs:
                    f_pdbs.write(response.text)
                    if args['--verbose']:
                        print("{} downloaded.".format(pdb_id.rstrip()))

    # Report of the download process
    print("All pdb files have been downloaded.")
    if nb_errors > 0:
        print("{} error(s) occurred.".format(nb_errors),
              "Watch the log_error.txt.")


def dssp(args):
    """
    Assign for each amino acid from protein in pdb file, a secondary structure.

    """
    # Check if the dssp' directory exists and create it if not
    if not os.path.isdir(DSSP_DIR):
        os.mkdir(DSSP_DIR)

    # Add secondary structures for each pdb file given in argument.
    for pdb in args['<pdb_file>']:
        nom_fichier = DSSP_DIR + os.path.splitext(os.path.basename(
            pdb))[0] + ".dssp"

        if args['--verbose']:
            subprocess.call([DSSP, '-vi', pdb, '-o', nom_fichier])
            print("Secondary struture assigned in {}.".format(nom_fichier))
        else:
            subprocess.call([DSSP, '-i', pdb, '-o', nom_fichier])


def convert(args):
    """Convert dssp file in csv and eleminate unnecessary informations."""
    if args['--rm']:
        try:
            os.remove(MUTUAL)
        except OSError:
            print("The file: {}, has already been removed or".format(MUTUAL),
                  "was never created.", file=sys.stderr)

    for dssp_file in args['<dssp_file>']:
        subprocess.call([DSSP_TO_CSV, dssp_file])

        if args['--verbose']:
            print("{} informations have been merged in mutual.csv.".format(
                os.path.basename(dssp_file)))


def get_angles_vectors():
    """
    Write a file of vectors (vectors.csv) composed of 8 angles per line
    with PHI and PSI angles extract from a previous file (angles.csv) generated
    by transforming chains of amino acids in chains of angles.
    Each chain is represented by on line in this file.
    """

    if not os.path.isfile("results/mutual.csv"):
        print("The mutual csv does not exists in results directory.",
              file=sys.stderr)
        print("Generate it by converting dssp files with convert command.",
              file=sys.stderr)

    # Creation of the angle.csv file
    # Reading the mutual.csv with specific options
    with open(MUTUAL, 'r') as f_mutual:
        reader = csv.reader(f_mutual, delimiter=';')

        # Opening the file that will contains PHI and PSI angle between each
        # pair of amino acids
        with open(ANGLES, 'w') as f_angles:
            # PHI -> field n°23  PSI -> field n°24
            for line in reader:
                # PHI and PSI == 360.0 --> chain break
                if line[23] == '360.0' and line[24] == '360.0':
                    pass
                # PHI == 360.0 --> beginnig of the polypeptide
                elif line[23] == '360.0':
                    f_angles.write(line[24] + ";")
                # PSI == 360.0 and PHI != 360.0 --> ending of the polypeptide
                elif line[24] == '360.0' and line[23] != '360.0':
                    f_angles.write(line[23] + "\n")
                else:
                    f_angles.write(line[23] + ";" + line[24] + ";")

    # Creation of the vectors.csv file by splitting chains of angles with
    # a sliding window of length 8.
    with open(VECTORS, 'w') as f_vectors:
        with open(ANGLES) as f_angles:
            for line in f_angles:
                list_angles = line.strip('\n').split(";")
                for i in range(0, (len(list_angles)-6), 2):
                    f_vectors.write(';'.join(list_angles[i:i+8])+"\n")


def main():
    """Launch the script"""
    # Parsing given arguments
    args = docopt.docopt(__doc__, version='prepare kohonen 0.0.1')

    if args['pdb']:
        if is_file(args['<ids_file>']):
            download_pdbs(args)
    elif args['dssp']:
        if verify_file(args['<pdb_file>'], '.pdb'):
            dssp(args)
    elif args['convert']:
        if verify_file(args['<dssp_file>'], '.dssp'):
            convert(args)
    elif args['angle']:
        get_angles_vectors()

if __name__ == '__main__':
    main()
