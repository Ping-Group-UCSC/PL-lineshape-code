#!/usr/bin/env python3

import os
import yaml
import re
import h5py
import numpy as np
import sys

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
        raise ValueError("Interface {} not implemented".format(interface))
    return wk, list_delta_r

def read_dynmat_mold_qe(nat, file):
    """
    Check if file exists and return freq and list_delta_r
    """
    if os.path.exists(file):
        with open(file) as f:
            print(indent, "Read dynmat mold from", file)
            lines = f.readlines()
    else:
        raise ValueError("couldn't find phonon file {}".format(file)) 
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
    interfaced with phonopy: read phonon info from band.yaml
    """
    #read data from file
    if not os.path.exists(f_band):
        raise ValueError("couldn't find phonon file {}".format(f_band))      
    if ".yaml" in f_band:
        wk, list_delta_r = read_phonon_yaml(f_band)
    elif ".hdf5" in f_band:
        wk, list_delta_r = read_phonon_hdf5(f_band)
    else: 
        raise ValueError("wrong format! need hdf5 or yaml phonon file.")
        
    return wk, list_delta_r

def read_phonon_yaml(f_band):
    #read phonon from phononpy file : bands.yaml
    with open(f_band, "r") as f:
        print("Read phonon from", f_band)
        ph_data = yaml.safe_load(f)  

    #read frequency and mode
    for d in ph_data['phonon']:
        if d['q-position']==[0.0, 0.0, 0.0]: 
            #the frequency in phonopy is THz; frequency of qe is cm-1;
            freq = [band['frequency']*THzToCm for band in d['band']]
            wk = [fre*conv_freq_to_omega for fre in freq]
            vec = [band['eigenvector'] for band in ph_data['phonon'][0]['band']] #phonon mode
            #we keep only the real part;
            list_delta_r = [
     [[r_component[0] for r_component in atoms] for atoms in mode ] for mode in vec
                            ]
            break
        raise ValueError("didn't find q=(0,0,0) point")
        
    return wk, list_delta_r

def read_phonon_hdf5(file):
    """
    Read the phonon freq. and eigenvector in hdf5 file from phononpy output
    file='band.hdf5'
    """
    #open the *.hdf5
    try :
        f = h5py.File(file,'r')
    except OSError:
        print("Could not open/read file: band.hdf5 or qpoints.hdf5")
        sys.exit()

    #check the file have content
    #print(list(f)) #--> ['dynamical_matrix', 'eigenvector', 'frequency', 'qpoint']
    if 'eigenvector' not in list(f):
        print("The hdf5 file doen't have 'eigenvector'!! re-run phonopy with '--eigvecs'")
        sys.exit()
    elif 'frequency' not in list(f):
        print("The hdf5 file doen't have 'frequency'!! Is it phonon data?")
        sys.exit()

    #read the content; eigenvector and frequeency
    tmp_eigenvector = f["eigenvector"][0].T #only Gamma point
    leng = len(tmp_eigenvector)
    tmp_eigenvector = tmp_eigenvector.reshape((leng, leng//3, 3)) #reshape the form; [freq.][atom][x,y,z]
    tmp_freq = f["frequency"][0].astype('float64') #only Gamma point
    #print("freq in band.yaml [THz]:")
    #print(frequencies)
    
    #sort the freq. and eigenv
    frequencies = np.sort(tmp_freq)
    frequencies_index = np.argsort(tmp_freq)
    del tmp_freq
    
    if list(frequencies_index) != range(leng):
        
        print("\n\tNow sort the freq. and eigenvector\n")
        eigenvector = [tmp_eigenvector[i] for i in frequencies_index]
        del tmp_eigenvector
        eigenvector = np.array(eigenvector)
    
    

    #when the freq is negative,
    print("\t freq inserted [THz] :",frequencies)
    count = np.where(frequencies <0)
    frequencies = np.where(frequencies >0, frequencies, 0 )
    print("\t # of negative freq :",len(count[0]))
    print()

    #change the format as code
    wk = (frequencies * THzToCm * conv_freq_to_omega) # array type & unit : rad/s
    list_delta_r = eigenvector # array type

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


def read_ZPL(f1, f2, interface="qe"):
    if os.path.exists(f1) and os.path.exists(f2):
        if interface=="qe":
            zpl = abs(read_totE(f1) - read_totE(f2))
            return zpl
        elif interface=="vasp":
            zpl = abs(read_totE_vasp(f1) - read_totE_vasp(f2))
            return zpl
    elif not os.path.exists(f1):
        print("The file %s does not exist" % f1)
        return None
    else:
        print("The file %s does not exist" % f2)
        return None
    
def read_totE_vasp(file):
    etot=None
    elec_conv=False
    ion_conv=False
    with open(file) as f:
        lines = f.readlines()
        for line in lines[::-1]:
            if "EDIFF is reached" in line:
                elec_conv= True
            if "reached required accuracy - stopping structural energy minimisation" in line:
                ion_conv = True
            m=re.search("TOTEN\s+=\s+(.+)\s+eV",line)
            if m:
                etot=m.groups()[-1]
                break
        print(f"electron convergence reached? {elec_conv}. Ionic convergence reached? {ion_conv} ")
    return float(etot)

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
