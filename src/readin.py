#!/bin/python3
'''
This script reads a Gaussian input file and extracts the basis set.
'''

# This script reads a Gaussian input file and extracts the basis set.
import periodictable as pt # type: ignore

def orcabasisformat(f):
    '''
    Reads an ORCA basis set file.
    '''
    # create a dictionary of the elements
    pesdict = {
        "symb": [],
        "numb": [],
        "startline": [],
        "endline": [],
        "nbf": [],
        "npr": [],
    }
    for el in pt.elements:
        if el.number == 0:
            continue
        elif el.number >= 87:
            break
        pesdict["symb"].append(el.symbol)
        pesdict["numb"].append(int(el.number))

    # read the file
    lines = f.readlines()
    for i in range(0, len(lines)):
        if "NewGTO" in lines[i].strip():
            pesdict["startline"].append(i + 1)
        if "end" in lines[i].strip():
            pesdict["endline"].append(i)

    bf = []
    exponents = []
    coefficients = []
    # create a dictionary with element symbol, number, nbf, npr, list[angular momenta],
    #   list[exponents], list[coefficients]
    basisset = {
        "symb": [],
        "numb": [],
        "nbf": [],
        "npr": [],
        "angmom": [],
        "lnpr": [],
        "exponents": [],
        "coefficients": [],
    }
    # go through all entries in the dictionary and read the basis functions
    for i in range(0, len(pesdict["symb"])):
        l = 0
        exptmp = []
        coefftmp = []
        bfdict = {"lnpr": [], "angmom": []}  # dictionary for the basis functions
        # set number of basis functions and primitives to zero for each element
        pesdict["nbf"].append(0)
        pesdict["npr"].append(0)
        for j in range(pesdict["startline"][i], pesdict["endline"][i]):
            if "S" in lines[j].strip():
                bfdict["angmom"].append("S")
            elif "P" in lines[j].strip():
                bfdict["angmom"].append("P")
            elif "D" in lines[j].strip():
                bfdict["angmom"].append("D")
            elif "F" in lines[j].strip():
                bfdict["angmom"].append("F")
            elif "G" in lines[j].strip():
                print("No support of G or higher functions.")
                exit()
            elif "end" in lines[j].strip():
                break
            else:
                continue
            num = int(lines[j].strip().split()[1])
            bfdict["lnpr"].append(num)
            pesdict["nbf"][i] = pesdict["nbf"][i] + 1
            pesdict["npr"][i] = pesdict["npr"][i] + bfdict["lnpr"][l]
            for m in range(j + 1, j + bfdict["lnpr"][l] + 1):
                exptmp.append(float(lines[m].strip().split()[1]))
                coefftmp.append(float(lines[m].strip().split()[2]))
            l += 1

        # Development printouts
        # print(pesdict["symb"][i],pesdict["numb"][i],pesdict["nbf"][i],pesdict["npr"][i])
        #######################
        exponents.append(exptmp)
        coefficients.append(coefftmp)
        bf.append(bfdict)

        # Development printouts
        # print(bf[i])
        # print(exponents[i])
        # print(coefficients[i])
        #######################

        # fill the dictionary with the basis set information
        basisset["symb"].append(pesdict["symb"][i])
        basisset["numb"].append(pesdict["numb"][i])
        basisset["nbf"].append(pesdict["nbf"][i])
        basisset["npr"].append(pesdict["npr"][i])
        basisset["angmom"].append(bf[i]["angmom"])
        basisset["lnpr"].append(bf[i]["lnpr"])
        basisset["exponents"].append(exptmp)
        basisset["coefficients"].append(coefftmp)

        # end of iteration

    return basisset


def orcaecpformat(f):
    '''
    Reads an ORCA ECP file.
    '''
    # create a dictionary of the elements
    startelemfound = False
    ecpdict = {
        "symb": [],
        "numb": [],
        "startline": [],
        "endline": [],
        "nbf": [],
        "npr": [],
        "ncore": [],
        "lmax": [],
    }
    lines = f.readlines()
    for el in pt.elements:
        if startelemfound is False:
            if el.symbol.lower() == lines[0].strip().split()[1].lower():
                startelemfound = True
            else:
                continue
        elif el.number >= 87:
            break
        ecpdict["symb"].append(el.symbol)
        ecpdict["numb"].append(int(el.number))

    # read the file
    for i in range(0, len(lines)):
        if "NewECP" in lines[i].strip():
            ecpdict["startline"].append(i + 1)
        if "end" in lines[i].strip():
            ecpdict["endline"].append(i)

    ecpf = []
    exponents = []
    coefficients = []
    ecpnfactor = []
    # create a dictionary with element symbol, number, nbf, npr, list[angular momenta],
    #   list[exponents], list[coefficients]
    ecpset = {
        "symb": [],
        "numb": [],
        "nbf": [],
        "npr": [],
        "angmom": [],
        "lnpr": [],
        "exponents": [],
        "coefficients": [],
        "ecpnfactor": [],
        "ncore": [],
        "lmax": [],
    }
    # go through all entries in the dictionary and read the basis functions
    for i in range(0, len(ecpdict["symb"])):
        l = 0
        exptmp = []
        coefftmp = []
        ecpnfactortmp = []
        bfdict = {"lnpr": [], "angmom": []}  # dictionary for the basis functions
        # set number of basis functions and primitives to zero for each element
        ecpdict["nbf"].append(0)
        ecpdict["npr"].append(0)
        ecpdict["ncore"].append(int(lines[ecpdict["startline"][i]].strip().split()[1]))
        ecpdict["lmax"].append(
            lines[ecpdict["startline"][i] + 1].strip().split()[1].lower()
        )
        for j in range(ecpdict["startline"][i] + 2, ecpdict["endline"][i]):
            if "s" in lines[j].strip():
                bfdict["angmom"].append("S")
            elif "p" in lines[j].strip():
                bfdict["angmom"].append("P")
            elif "d" in lines[j].strip():
                bfdict["angmom"].append("D")
            elif "f" in lines[j].strip():
                bfdict["angmom"].append("F")
            elif "g" in lines[j].strip():
                bfdict["angmom"].append("G")
            elif "h" in lines[j].strip():
                bfdict["angmom"].append("H")
            elif "i" in lines[j].strip():
                print("No support of i or higher functions.")
                exit()
            elif "end" in lines[j].strip():
                break
            else:
                continue
            num = int(lines[j].strip().split()[1])
            bfdict["lnpr"].append(num)
            ecpdict["nbf"][i] = ecpdict["nbf"][i] + 1
            ecpdict["npr"][i] = ecpdict["npr"][i] + bfdict["lnpr"][l]
            for m in range(j + 1, j + bfdict["lnpr"][l] + 1):
                exptmp.append(float(lines[m].strip().split()[1]))
                coefftmp.append(float(lines[m].strip().split()[2]))
                ecpnfactortmp.append(int(lines[m].strip().split()[3]))
            l += 1

        # Development printouts
        # print(pesdict["symb"][i],pesdict["numb"][i],pesdict["nbf"][i],pesdict["npr"][i])
        #######################
        exponents.append(exptmp)
        coefficients.append(coefftmp)
        ecpnfactor.append(ecpnfactortmp)
        ecpf.append(bfdict)

        # Development printouts
        # print(bf[i])
        # print(exponents[i])
        # print(coefficients[i])
        #######################

        # fill the dictionary with the basis set information
        ecpset["symb"].append(ecpdict["symb"][i])
        ecpset["numb"].append(ecpdict["numb"][i])
        ecpset["nbf"].append(ecpdict["nbf"][i])
        ecpset["npr"].append(ecpdict["npr"][i])
        ecpset["angmom"].append(ecpf[i]["angmom"])
        ecpset["lnpr"].append(ecpf[i]["lnpr"])
        ecpset["exponents"].append(exptmp)
        ecpset["coefficients"].append(coefftmp)
        ecpset["ecpnfactor"].append(ecpnfactortmp)
        ecpset["ncore"].append(ecpdict["ncore"][i])
        ecpset["lmax"].append(ecpdict["lmax"][i])

        # end of iteration

    return ecpset
