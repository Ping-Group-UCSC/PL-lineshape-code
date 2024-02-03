#!/usr/bin/env python3

import os

from constant import indent, Ry2eV, conv_freq_to_omega, THzToCm

def read_dynmat_mold(nat, file, interface='qe'):
    """
    interface: 'qe' , or 'phonopy'
    """
    if interface=='qe':
        wk, list_delta_r = read_dynmat_mold_qe(nat, file)
    elif interface=='phonopy':
        wk, list_delta_r = read_dynmat_mold_phonopy(file)
    else:
        raise ValueError("Interface not implemented")
    return wk, list_delta_r

def read_dynmat_mold_qe(nat, file):
    """
    Check if file exists and return freq and list_delta_r
    """
    if os.path.exists(file):
        with open(file) as f:
            print(indent, "Read dynmat mold from", file)
            lines = f.readlines()
    #    nat, nmodes = calc_nat_nmodes(lines)
    nmodes = nat * 3
    for i, line in enumerate(lines):
        if "[FREQ]" in line:
            # save next nmodes line to the list wk
            wk = [
                float(line.split()[0]) * conv_freq_to_omega
                for line in lines[i + 1 : i + 1 + nmodes]
            ]
        elif "[FR-NORM-COORD]" in line:
            list_delta_r = []
            begin = i + 2
            for count in range(nmodes):
                end = begin + nat
                # each part is a nat rows by 3 col array
                part = [
                    [float(val) for val in line.split()] for line in lines[begin:end]
                ]
                # list_delta_r stores all parts in a list of length nmodes
                list_delta_r.append(part)
                begin = end + 1
    return wk, list_delta_r

def read_dynmat_mold_phonopy(f_band):
    """
    Check if file exists and return freq and list_delta_r
    f_band:  band.yaml include band information;
    interfaced with phonopy
    """
    #read data from file
    if os.path.exists(f_band):
        import yaml
        with open(f_band, "r") as f:
            print("Read phonon from", f_band)
            data = yaml.safe_load(f)

    #read frequency and mode
    for d in data['phonon']:
        if d['q-position']==[0.0, 0.0, 0.0]: 
            #the frequency in phonopy is THz; frequency of qe is cm-1;
            freq = [band['frequency']*THzToCm for band in d['band']]
            wk = [fre*conv_freq_to_omega for fre in freq]
            vec = [band['eigenvector'] for band in data['phonon'][0]['band']] #phonon mode
            #we keep only the real part;
            list_delta_r = [
     [[r_component[0] for r_component in atoms] for atoms in mode ] for mode in vec
                            ]
            break
        raise ValueError("didn't find q=(0,0,0) point")
        
    return wk, list_delta_r

def calc_nat_nmodes(lines):
    count_flag = False
    for line in lines:
        if "[FREQ]" in line:
            count = 0
            count_flag = True
            continue
        elif "[FR-COORD]" in line:
            break
        elif count_flag:
            count += 1
    nmodes = int(count)
    nat = int(nmodes / 3)
    return nat, nmodes


def read_ZPL(f1, f2):
    if os.path.exists(f1) and os.path.exists(f2):
        zpl = abs(read_totE(f1) - read_totE(f2))
        return zpl
    elif not os.path.exists(f1):
        print("The file %s does not exist" % f1)
        return None
    else:
        print("The file %s does not exist" % f2)
        return None


def read_totE(file):
    with open(file) as f:
        lines = f.readlines()
        for line in lines[::-1]:
            if "total energy" in line and "is the sum" not in line:
                tag = line[:2].strip()
                break
        if tag != "!" and tag != "!!":
            print("Cannot recognize total energy in %s" % file.replace(".in", ".out"))
        for line in lines[::-1]:
            if line.startswith(tag):
                etot = float(line.split()[-2]) * Ry2eV
                break
    return etot
