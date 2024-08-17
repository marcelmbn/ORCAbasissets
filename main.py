#!/usr/bin/env python3
'''
This program processes basis sets and ECPs.
'''

import argparse
from src.readin import orcabasisformat, orcaecpformat
from src.writebasis import orcabasissetcode, orcaecpcode, xtb_tblite_format_basis, turbomole_format_basis, turbomole_format_ecp


def printbasis(basis, desiredelem):
    '''
    Print the basis set to the screen in a proprietary format.
    '''
    # print only the full basis function for the desired element
    print("The number basis functions for the desired element is:")
    print(basis["nbf"][desiredelem - 1])
    print("The number of primitives for the desired element is:")
    print(basis["npr"][desiredelem - 1])
    print("The basis functions are:")
    k = 0
    for i in range(0, basis["nbf"][desiredelem - 1]):
        print("Basis function " + str(i + 1) + ":")
        print("  Angular momentum: " + basis["angmom"][desiredelem - 1][i])
        print("  Number of primitives: " + str(basis["lnpr"][desiredelem - 1][i]))
        print("  Exponents and coefficients:")
        for j in range(0, basis["lnpr"][desiredelem - 1][i]):
            print(j + 1)
            print("    " + str(basis["exponents"][desiredelem - 1][k]))
            print("    " + str(basis["coefficients"][desiredelem - 1][k]))
            k += 1


def printecp(ecp, desiredelem):
    '''
    Print the ECP to the screen in a proprietary format.
    '''
    # print only the full basis function for the desired element
    print(ecp["numb"][0])
    if ecp["numb"][0] > desiredelem:
        print("The element of interest is not treated with ECP.")
        exit()
    desiredelem = desiredelem - ecp["numb"][0]  # subtract elements treated without ECP
    print("Element of interest:")
    print("  Symbol: " + ecp["symb"][desiredelem])
    print(
        "Maximum ang. momentum and n_core "
        + str(ecp["lmax"][desiredelem])
        + " "
        + str(ecp["ncore"][desiredelem])
    )
    print("The number of ECP functions for the desired element is:")
    print(ecp["nbf"][desiredelem])
    print("The number of primitives for the desired element is:")
    print(ecp["npr"][desiredelem])
    print("The ECP functions are:")
    k = 0
    for i in range(0, ecp["nbf"][desiredelem]):
        print("ECP function " + str(i + 1) + ":")
        print("  Angular momentum: " + ecp["angmom"][desiredelem][i])
        print("  Number of primitives: " + str(ecp["lnpr"][desiredelem][i]))
        print("  Exponents and coefficients:")
        for j in range(0, ecp["lnpr"][desiredelem][i]):
            print(f"  {j+1} {ecp['exponents'][desiredelem][k]} {ecp['coefficients'][desiredelem][k]} {ecp['ecpnfactor'][desiredelem][k]}")
            k += 1

# check for command line arguments based on argparse
parser = argparse.ArgumentParser(
    description="Read in a Gaussian input file and extract the basis set."
)
parser.add_argument(
    "--ifile", metavar="inputfile", type=str, nargs=1, help="the input file"
)
parser.add_argument(
    "--element",
    metavar="printelement",
    type=int,
    nargs=1,
    help="the element for which the basis functions should be printed",
)
# check optionally if ECP or basis set mode using the argument --ecp and --basis
parser.add_argument(
    "--ecpmode",
    action="store_true",
    help="use this flag if the input file contains an ECP",
    default=False,
)
parser.add_argument(
    "--basismode",
    action="store_true",
    help="use this flag if the input file contains a basis set",
    default=False,
)
parser.add_argument(
    "--format",
    type=str,
    help="the format of the basis set",
    default="orca",
    required=False,
)
parser.add_argument(
    "-v", "--verbose", action="store_true", help="use this flag for extreme output"
)

args = parser.parse_args()

# check if the input argument is given
try:
    str(args.ifile[0])
except:
    print("No input argument.")
    exit()

# check if the input file can be opened
try:
    file = open(str(args.ifile[0]), "r", encoding="utf-8")
    print("File " + str(args.ifile[0]) + " opened.")
except:
    print("Could not open file " + str(args.ifile[0]) + ".")
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
    ecpmode = False
    basismode = True
elif args.basismode and not args.ecpmode:
    print("Basis set mode.")
    ecpmode = False
    basismode = True
elif not args.basismode and args.ecpmode:
    print("ECP mode.")
    ecpmode = True
    basismode = False
else:
    print("Both modes given.\nAssuming basis set mode.")
    ecpmode = False
    basismode = True

if basismode:
    # read the file
    # basis is a dictionary with the following keys and values:
    # symb: list of the symbols of the elements
    # numb: list of the atomic numbers of the element
    basis = orcabasisformat(file)
    if args.verbose:
        printbasis(basis, desiredelem)

    print("The format of the basis set is:")
    print(args.format)
    if args.format == "orca":
        orcabasissetcode(basis)
    elif args.format == "xtb":
        xtb_tblite_format_basis(basis)
    elif args.format == "turbomole":
        turbomole_format_basis(basis)
    else:
        print("Format not supported.")
        exit()

elif ecpmode:
    # read the file
    # ecp is a dictionary with the following keys and values:
    # symb: list of the symbols of the elements
    # numb: list of the atomic numbers of the element
    ecp = orcaecpformat(file)
    if args.verbose:
        printecp(ecp, desiredelem)

    # if args.format not orca
    if args.format == "orca":
        orcaecpcode(ecp)
    if args.format == "turbomole":
        turbomole_format_ecp(ecp)
    else:
        print("Format not supported.")
        exit()
else:
    print("Something went wrong.")
    exit()
