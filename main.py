#!/bin/python3

import argparse
from src.readin import orcaformat
from src.writebasis import orcabasissetcode

def printbasis(basis,desiredelem):
    # print only the full basis function for the desired element
    print("The number basis functions for the desired element is:")
    print(basis["nbf"][desiredelem-1])
    print("The number of primitives for the desired element is:")
    print(basis["npr"][desiredelem-1])
    print("The basis functions are:")
    k=0
    for i in range(0,basis["nbf"][desiredelem-1]):
        print("Basis function "+str(i+1)+":")
        print("  Angular momentum: "+basis["angmom"][desiredelem-1][i])
        print("  Number of primitives: "+str(basis["lnpr"][desiredelem-1][i]))
        print("  Exponents and coefficients:")
        for j in range(0,basis["lnpr"][desiredelem-1][i]):
            print(j+1)
            print("    "+str(basis["exponents"][desiredelem-1][k]))
            print("    "+str(basis["coefficients"][desiredelem-1][k]))
            k+=1

basismode = False
ecpmode = False

# check for command line arguments based on argparse
parser = argparse.ArgumentParser(description='Read in a Gaussian input file and extract the basis set.')
parser.add_argument("--ifile", metavar='inputfile', type=str, nargs=1, help='the input file')
parser.add_argument("--element", metavar='printelement', type=int, nargs=1, help='the element for which the basis functions should be printed')
# check optionally if ECP or basis set mode using the argument --ecp and --basis
parser.add_argument("--ecpmode", action='store_true', help='use this flag if the input file contains an ECP')
parser.add_argument("--basismode", action='store_true', help='use this flag if the input file contains a basis set')
parser.add_argument("-v","--verbose", action='store_true', help='use this flag for extreme output')

args = parser.parse_args()

# check if the input argument is given
try:
    str(args.ifile[0])
except:
    print("No input argument.")
    exit()

# check if the input file can be opened
try:
    file = open(str(args.ifile[0]))
    print("File "+str(args.ifile[0])+" opened.")
except:
    print("Could not open file "+str(args.ifile[0])+".")
    exit()

# check if the input argument is given
try:
    str(args.element[0])
    desiredelem = args.element[0]
except:
    print("No desired element to print.")
    desiredelem = 86
    
# check if basisset or ecp mode is given

if not args.basismode and not args.ecpmode:
    print("No mode given.\nAssuming basis set mode.")
    basismode = True
elif args.basismode and not args.ecpmode:
    print("Basis set mode.")
    basismode = True
elif not args.basismode and args.ecpmode:
    print("ECP mode.")
    ecpmode = True
else:
    print("Both modes given.\nAssuming basis set mode.")
    basismode = True

if basismode:
    # read the file
    # basis is a dictionary with the following keys and values:
    # symb: list of the symbols of the elements
    # numb: list of the atomic numbers of the element
    basis = orcaformat(file)
    if args.verbose:
        printbasis(basis,desiredelem)

    orcabasissetcode(basis)
