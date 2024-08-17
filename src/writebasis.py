"""
This module contains functions to write the basis set in the format of
different quantum chemistry programs.
"""

import os
import numpy as np

# dictionary for relationship between angular momenta
# in letters ("S","P","D","F") and numbers (0,1,2,3)
angmomdict = {"S": 0, "P": 1, "D": 2, "F": 3, "G": 4, "H": 5, "I": 6}


def orcabasissetcode(bas):
    """
    Write the basis set in the format of the ORCA input file
    """
    # Write the basis set in the format of the ORCA input file
    # The input is a dictionary with the following keys and values:
    path = "output"
    # Check whether the specified path exists or not
    isExist = os.path.exists(path)
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs(path)
        print("Output directory was created!")

    ofile = open("output/basis_orcasource.txt", "w", encoding="utf-8")
    for i in range(0, len(bas["symb"])):
        l = 0
        ofile.write("  // --------------------------------------------------------\n")
        ofile.write(
            "  // element no " + str(bas["numb"][i]) + " (" + bas["symb"][i] + ")\n"
        )
        ofile.write("  // --------------------------------------------------------\n")
        ofile.write(
            "  B.Init({0:4d},{1:4d});\n".format(
                int(bas["numb"][i]), int((bas["nbf"][i]))
            )
        )
        for j in range(0, bas["nbf"][i]):
            ofile.write("  // Basis function " + bas["angmom"][i][j] + "\n")
            ofile.write(
                "  BG[{0:3d}][{1:3d}].l     ={2:2d};\n".format(
                    bas["numb"][i], j, angmomdict[bas["angmom"][i][j]]
                )
            )
            ofile.write(
                "  BG[{0:3d}][{1:3d}].ng    ={2:2d};\n".format(
                    bas["numb"][i], j, bas["lnpr"][i][j]
                )
            )
            for k in range(0, bas["lnpr"][i][j]):
                ofile.write(
                    "  BG[{0:3d}][{1:3d}].a[{2:2d}] ={3:25.10f};   \
BG[{0:3d}][{1:3d}].d[{2:2d}] ={4:25.10f};\n".format(
                        bas["numb"][i],
                        j,
                        k,
                        bas["exponents"][i][l],
                        bas["coefficients"][i][l],
                    )
                )
                l += 1
        ofile.write("\n\n")
    ofile.close()


def orcaecpcode(ecp):
    """
    Write the ECP in the format of the ORCA input file
    """
    # Write the ECP in the format of the ORCA input file
    # The input is a dictionary with the following keys and values:
    path = "output"
    # Check whether the specified path exists or not
    isExist = os.path.exists(path)
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs(path)
        print("Output directory was created!")

    ofile = open("output/ecp_orcasource.txt", "w", encoding="utf-8")
    for i in range(0, len(ecp["symb"])):
        l = 0
        ofile.write("  // --------------------------------------------------------\n")
        ofile.write(
            "  // element no " + str(ecp["numb"][i]) + " (" + ecp["symb"][i] + ")\n"
        )
        ofile.write("  // --------------------------------------------------------\n")
        ofile.write(
            "  E->U[ {0:3d} ] = new ecpPotential;\n".format(int(ecp["numb"][i]))
        )
        ofile.write(
            "  E->U[ {0:3d} ]->initialize({1:4d},{2:4d},{3:4d} );\n".format(
                ecp["numb"][i],
                int(angmomdict[ecp["lmax"][i].upper()]),
                ecp["ncore"][i],
                ecp["nbf"][i],
            )
        )
        ofile.write("  E->U[ {0:3d} ]->element    ={0:4d};\n".format(ecp["numb"][i]))
        ofile.write(
            "  CopyECPCommentary( E->U[ {0:3d} ], \
all_ecpname, all_comment, all_citation, all_source );\n".format(
                ecp["numb"][i]
            )
        )
        ofile.write("  E->U[ {0:3d} ]->cpp        = NULL;\n".format(ecp["numb"][i]))
        ofile.write("  EG = E->U[ {0:3d} ]->U_l;\n\n".format(ecp["numb"][i]))
        for j in range(0, ecp["nbf"][i]):
            ofile.write("  EG[{0:4d} ].l ={0:2d};\n".format(j))
            ofile.write("  EG[{0:4d} ].K ={1:2d};\n".format(j, ecp["lnpr"][i][j]))
            for k in range(0, ecp["lnpr"][i][j]):
                ofile.write(
                    "  EG[{0:4d} ].p[{1:4d} ].a ={2:23.12f} ;".format(
                        j, k, ecp["exponents"][i][l]
                    )
                )
                ofile.write(
                    "  EG[{0:4d} ].p[{1:4d}  ].d ={2:23.12f} ;".format(
                        j, k, ecp["coefficients"][i][l]
                    )
                )
                ofile.write(
                    "  EG[{0:4d} ].p[{1:4d}  ].n ={2:4.1f} ;\n".format(
                        j, k, float(ecp["ecpnfactor"][i][l])
                    )
                )
                l += 1
        ofile.write("\n\n")


def xtb_tblite_format_basis(bas):
    """
    Write the basis set in the format of the xtb tblite file
    """
    path = "output"
    isExist = os.path.exists(path)
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs(path)
        print("Output directory was created!")

    # find out the maximal number of primitive functions per basis function
    # in the whole basis set
    maxnpr = 0
    for i in range(0, len(bas["nbf"])):
        for j in range(0, bas["nbf"][i]):
            if bas["lnpr"][i][j] > maxnpr:
                maxnpr = bas["lnpr"][i][j]
    print("Maximum number of primitive functions per basis function: " + str(maxnpr))
    maxshell = 7

    # open("output/basis_xtbsource.txt", "w", encoding="utf-8")
    with open("output/basis_xtbsource.txt", "w", encoding="utf-8") as ofile:
        # iterate over all elements until Z=86
        print("----------- EXPONENTS -----------")
        for i in range(0, len(bas["symb"])):
            if bas["numb"][i] > 57 and bas["numb"][i] < 72:
                continue
            print("Element no " + str(bas["numb"][i]) + " (" + bas["symb"][i] + ")")
            # write the basis set in Fortran format as follows:
            # create one Fortran array with the exponents
            # and one with the coefficients for each element
            # the arrays should have max(nprim) x max(nbf) elements
            # and are setup with zeros for the unused elements
            # and using the Fortran "reshape" function to create the
            # correct dimensions
            exponents = np.zeros((maxnpr, maxshell), dtype=float)
            # go through all exponents of the basis functions of the element i
            m = 0
            for j in range(0, bas["nbf"][i]):
                for k in range(0, bas["lnpr"][i][j]):
                    exponents[k][j] = bas["exponents"][i][m]
                    m += 1

            # print the numpy array in nice format to screen
            print(exponents)

            # write the exponents in Fortran format
            print(f"exponents(:, :, {i+1:2d}) = reshape([&", file=ofile)
            print(exponents.shape[0], exponents.shape[1])
            for j in range(exponents.shape[1]):
                print("& ", file=ofile, end="")
                for k in range(exponents.shape[0] - 1):
                    print(f"{exponents[k][j]:15.10f}_wp, ", file=ofile, end="")
                if j == exponents.shape[1] - 1:
                    print(
                        f"{exponents[-1][j]:15.10f}_wp], (/max_prim, max_shell/))",
                        file=ofile,
                    )
                else:
                    print(f"{exponents[-1][j]:15.10f}_wp, &", file=ofile)
            print("", file=ofile)

        print("----------- COEFFICIENTS -----------")
        for i in range(0, len(bas["symb"])):
            if bas["numb"][i] > 57 and bas["numb"][i] < 72:
                continue
            print("Element no " + str(bas["numb"][i]) + " (" + bas["symb"][i] + ")")
            coefficients = np.zeros((maxnpr, maxshell), dtype=float)
            # go through all coefficients of the basis functions of the element i
            m = 0
            for j in range(0, bas["nbf"][i]):
                print("Basis function " + str(j) + ":")
                for k in range(0, bas["lnpr"][i][j]):
                    coefficients[k][j] = bas["coefficients"][i][m]
                    m += 1

            # print the numpy array in nice format to screen
            print(coefficients)

            # write the coefficients in Fortran format
            print(f"coefficients(:, :, {i+1:2d}) = reshape([&", file=ofile)
            for j in range(coefficients.shape[1]):
                print("& ", file=ofile, end="")
                for k in range(coefficients.shape[0] - 1):
                    print(f"{coefficients[k][j]:15.10f}_wp, ", file=ofile, end="")
                if j == coefficients.shape[1] - 1:
                    print(
                        f"{coefficients[-1][j]:15.10f}_wp], (/max_prim, max_shell/))",
                        file=ofile,
                    )
                else:
                    print(f"{coefficients[-1][j]:15.10f}_wp, &", file=ofile)
            print("", file=ofile)

        # print the number of basis functions and primitives for each element
        print("----------- NUMBER OF BASIS FUNCTIONS -----------", file=ofile)
        print("\ninteger, parameter :: nshell(max_elem) = [ &", file=ofile)
        l = 0
        for i in range(0, len(bas["symb"])):
            print("Element no " + str(bas["numb"][i]) + " (" + bas["symb"][i] + ")")
            print("The number basis functions for the desired element is:")
            print(bas["nbf"][i])
            if bas["nbf"][i] > 7:
                tmpbasnbf = 0
            else:
                tmpbasnbf = bas["nbf"][i]
            l += 1
            if l >= 20:
                print(f"{tmpbasnbf}, &", file=ofile)
                l = 0
                continue
            elif i >= len(bas["symb"]) - 1:
                print(f"{tmpbasnbf}]", file=ofile)
            elif l <= 1:
                print(f"& {tmpbasnbf}, ", file=ofile, end="")
            else:
                print(f"{tmpbasnbf}, ", file=ofile, end="")

        # print the number of primitives for each basis function of each element
        print("----------- NUMBER OF PRIMITIVES -----------", file=ofile)
        print(
            "\ninteger, parameter :: n_prim(highest_elem, max_shell) = reshape([&",
            file=ofile,
        )
        l = 0
        for i in range(0, len(bas["symb"])):
            print("Element no " + str(bas["numb"][i]) + " (" + bas["symb"][i] + ")")
            print("The number of primitives for the desired element is:")
            print(bas["npr"][i])
            for j in range(0, maxshell):
                if bas["nbf"][i] > 7:
                    tmpbasnpr = 0
                else:
                    if j >= bas["nbf"][i]:
                        tmpbasnpr = 0
                    else:
                        tmpbasnpr = bas["lnpr"][i][j]
                l += 1
                if l >= 21:
                    print(f"{tmpbasnpr}, & ! up to element: {i+1}", file=ofile)
                    l = 0
                elif i >= len(bas["symb"]) - 1 and j >= maxshell - 1:
                    print(f"{tmpbasnpr}]", file=ofile)
                elif l <= 1:
                    print(f"& {tmpbasnpr}, ", file=ofile, end="")
                else:
                    print(f"{tmpbasnpr}, ", file=ofile, end="")

        # print the angular momentum for each basis function of each element
        print("----------- ANGULAR MOMENTUM -----------", file=ofile)
        print(
            "\ninteger, parameter :: angmom(max_elem, max_shell) = reshape([&",
            file=ofile,
        )
        l = 0
        for i in range(0, len(bas["symb"])):
            print("Element no " + str(bas["numb"][i]) + " (" + bas["symb"][i] + ")")
            print("The angular momentum for the desired element is:")
            print(bas["angmom"][i])
            for j in range(0, maxshell):
                if bas["nbf"][i] > 7:
                    tmpbasangmom = 0
                else:
                    if j >= bas["nbf"][i]:
                        tmpbasangmom = 0
                    else:
                        tmpbasangmom = angmomdict[bas["angmom"][i][j]]
                l += 1
                if l >= 21:
                    print(f"{tmpbasangmom}, & ! up to element: {i+1}", file=ofile)
                    l = 0
                elif i >= len(bas["symb"]) - 1 and j >= maxshell - 1:
                    print(f"{tmpbasangmom}]", file=ofile)
                elif l <= 1:
                    print(f"& {tmpbasangmom}, ", file=ofile, end="")
                else:
                    print(f"{tmpbasangmom}, ", file=ofile, end="")


# write the basis set in the following TURBOMOLE format:
# $basis
# *
# rb def2-TZVP
# *
#     2   s
#       7.4744618040           0.26997866363
#       6.7296180594          -0.42629251814
#     1   s
#       2.7816640004           1.0000000
#     1   s
#       0.53452175148          1.0000000
#     1   s
#       0.22368793034          1.0000000
#     1   s
#       0.32410407052D-01      1.0000000
#     1   s
#       0.14171047424D-01      1.0000000
#     4   p
#       5.6720643194           0.48114224135D-01
#       3.3320183956          -0.18485131426
#       0.80150054910          0.42811864954
#       0.36302220227          0.58673165411
#     1   p
#       0.15733924392          1.0000000
#     1   p
#       0.40000000000D-01      1.0000000
#     1   p
#       0.16000000000D-01      1.0000000
#     1   d
#       0.25907866956          0.76806746340D-01
#     1   d
#       0.42507438045D-01      0.37846160487
#     1   d
#       0.11909276840D-01      1.0000000
# *
# sr def2-TZVP
# *
#     2   s
#      10.000000000           -0.18530550262
#       8.5000000000           0.33970355376
#     1   s
#       3.0057048856           1.0000000
#     1   s
#       0.61161287650          1.0000000
#     1   s
#       0.27393841217          1.0000000
#     1   s
#       0.57435564563D-01      1.0000000
#     1   s
#       0.23338198665D-01      1.0000000
#     4   p
#       7.5883077869           0.33731690287D-01
#       3.6731307392          -0.20523185005
#       0.90496618455          0.49209972665
#       0.43310256408          0.62105296512
#     1   p
#       0.20222168964          1.0000000
#     1   p
#       0.72000000000D-01      1.0000000
#     1   p
#       0.25000000000D-01      1.0000000
#     3   d
#       3.6180810000          -0.75010000000D-02
#       0.99665600000          0.10809800000
#       0.39073500000          0.27854000000
#     1   d
#       0.12277000000          1.0000000
#     1   d
#       0.36655000000D-01      1.0000000
# *
# $end

def turbomole_format_basis(bas) :
    """
    Write the basis set in the format of the TURBOMOLE input file
    """
    path = "output"
    isExist = os.path.exists(path)
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs(path)
        print("Output directory was created!")

    ofile = open("output/basis_turbomolesource.txt", "w", encoding="utf-8")
    ofile.write("$basis\n*\n")
    for i in range(0, len(bas["symb"])):
        l = 0
        ofile.write(f"{bas['symb'][i].lower()} vDZP\n*\n")
        for j in range(0, bas["nbf"][i]):
            ofile.write("    " + str(bas["lnpr"][i][j]) + "   " + bas["angmom"][i][j].lower() + "\n")
            for _ in range(0, bas["lnpr"][i][j]):
                ofile.write(f"{bas["exponents"][i][l]:>25}{bas["coefficients"][i][l]:>25}\n")
                l += 1
        ofile.write("*\n")
    ofile.write("$end\n")
    ofile.close()

# writh the ECP in the following TURBOMOLE format:
# $ecp
# *
# rb vDZP-ecp
# *
#   ncore = 28   lmax = 3
# f
#      -12.3169000      2       3.8431140
# s-f
#       89.5001980      2       5.0365510
#        0.4937610      2       1.9708490
#       12.3169000      2       3.8431140
# p-f
#       58.5689740      2       4.2583410
#        0.4317910      2       1.4707090
#       12.3169000      2       3.8431140
# d-f
#       26.2248980      2       3.0231270
#        0.9628390      2       0.6503830
#       12.3169000      2       3.8431140
# *
# sr vDZP-ecp
# *
#   ncore = 28   lmax = 3
# f
#      -15.8059920      2       4.6339750
# s-f
#      135.4794300      2       7.4000740
#       17.5344630      2       3.6063790
#       15.8059920      2       4.6339750
# p-f
#       88.3597090      2       6.4848680
#       15.3943720      2       3.2880530
#       15.8059920      2       4.6339750
# d-f
#       29.8889870      2       4.6228410
#        6.6594140      2       2.2469040
#       15.8059920      2       4.6339750
# *
# $end

def turbomole_format_ecp(ecp):
    """
    Write the ECP in the format of the TURBOMOLE input file
    """
    path = "output"
    isExist = os.path.exists(path)
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs(path)
        print("Output directory was created!")

    ofile = open("output/ecp_turbomolesource.txt", "w", encoding="utf-8")
    ofile.write("$ecp\n*\n")
    for i in range(0, len(ecp["symb"])):
        l = 0
        ofile.write(f"{ecp['symb'][i].lower()} vDZP-ecp\n*\n")
        ofile.write(f"  ncore = {ecp['ncore'][i]:2d}   lmax = {angmomdict[ecp['lmax'][i].upper()]:2d}\n")
        ofile.write(f"{ecp['angmom'][i][-1].lower()}\n")
        for j in range(0, ecp["nbf"][i]):
            for _ in range(0, ecp["lnpr"][i][j]):
                if j == ecp["nbf"][i]-1:
                    ofile.write(f"{ecp['coefficients'][i][l]:>15}{ecp['ecpnfactor'][i][l]:>5}{ecp['exponents'][i][l]:>15}\n")
                l += 1
        l = 0
        for j in range(0, ecp["nbf"][i]-1):
            ofile.write(f"{ecp['angmom'][i][j].lower()}-{ecp['angmom'][i][ecp['nbf'][i]-1].lower()}\n")
            for _ in range(0, ecp["lnpr"][i][j]):
                ofile.write(f"{ecp['coefficients'][i][l]:>15}{ecp['ecpnfactor'][i][l]:>5}{ecp['exponents'][i][l]:>15}\n")
                l += 1
        ofile.write("*\n")
    ofile.write("$end\n")
    ofile.close()
