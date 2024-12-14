#!/usr/bin/env python3

import os
import numpy as np

from constant import Electron2Coulomb, indent, h_eVs


class inp_class:
    """
    class for organizing input
    """

    """
    defaults
    """
    # path to qe output: relax-gs, relax-cdftup1, relax-gs/PH
    path_gs = "./relax-gs/relax"
    path_ex = "./relax-cdftup1/relax"
    # phonon interface = qe or phonopy;
    phonon_interface = "qe"
    # phonon file, default is from qe:
    file_phonon = os.path.join(".", "relax-gs/PH/dynmat.mold")
    # zero-phonon line, units eV (converted to Joules)
    zpl = None
    # file to read previously calculated wk, Sk
    skfile = None
    # gamma factor in PL, units of s^-1 [Hz]
    gamma = 0.28571428571e14
    # smearing for S(hw), units eV (converted to Joules)
    smear = 0.006 * Electron2Coulomb
    # min energy for PL spectrum ([eV] --> converted [J])
    hw_min = 1.0 * Electron2Coulomb
    # max energy for PL spectrum ([eV] --> converted [J])
    hw_max = 2.1 * Electron2Coulomb
    # number of points to produce
    hw_steps = None
    resolution = 1000
    #integral options:
    cal_method = "FFT" #method for compute PL spectrum A_hw; 'integral' or 'FFT'
    #if cal_method = 'integral', then following options choose integral method and threshold
    integral_method = "quad_vec"    #integral method for calculating A_hw ('quad_vec' & 'romberg')
    # improves accuracy of integral, unitless
    tolerance = 1e-18
    # limit used in integral, units of s
    limit = 3.5e-13
    #options : hard coded, not read by pl.in

    #plot option
    PL_hw_xmin = 1.2 #[eV]
    PL_hw_xmax = 2.4 #[eV]
    #S_hw plot option
    S_hw_xmin = 0 #[eV]
    S_hw_xmax = 0.2 #[eV]
    plot_flag = True    #option for ploting all data
    plot_pHR = True     #option for ploting HR factor
    parallel = False    #option for integral with parallel (only for 'quad_vec')
    smear_end = None
    

    def gen_hw_array(self):
        if self.hw_steps != None:
            dhw = (self.hw_max - self.hw_min) / self.hw_steps # (hw_max - hw_min)/hw_step
        else:
            dhw = (Electron2Coulomb) /self.resolution # (1eV-0eV)/resolution

        self.hw_array = np.arange(self.hw_min, self.hw_max, dhw) #converted to [J]

    def contents(self):
        # return tuple of all inp contents
        return (
            self.path_gs,
            self.path_ex,
            self.phonon_interface,
            self.file_phonon,
            self.zpl,
            self.skfile,
            self.smear,
            self.limit,
            self.gamma,
            self.tolerance,
            self.gen_hw_array(),
            self.cal_method,
            self.integral_method
        )

    def print_contents(self):
        E2C = Electron2Coulomb
        print( "{} path_gs   = {}".format(2 * indent, self.path_gs) )
        print( "{} path_ex   = {}".format(2 * indent, self.path_ex) )
        print( "{} phonon_interface   = {}".format(2 * indent, self.phonon_interface) )
        print( "{} file_phonon   = {}".format(2 * indent, self.file_phonon))
        if 'zpl' in vars():
             print( "{} zpl (eV)     = {}".format(2 * indent, str(self.zpl / E2C)) )
        else:
             print( "{} zpl (eV)     = {}".format(2 * indent, "from input") )
        print( "{} skfile       = {}".format(2 * indent, self.skfile) )
        print( "{} smear (eV)   = {:10.6e}".format(2 * indent, self.smear / E2C) )
        print( "{} limit (s)    = {:10.6e}".format(2 * indent, self.limit) )
        print( "{} gamma (s^-1) = {:10.6e}".format(2 * indent, self.gamma) )
        print( "{} gamma (eV) = {:10.6e} (E=hf)".format(2 * indent, self.gamma*h_eVs) )
        print( "{} tolerance    = {:10.6e}".format(2 * indent, self.tolerance) )
        print( "{} hw_min (eV)  = {:10.6f}".format(2 * indent, self.hw_min / E2C) )
        print( "{} hw_max (eV)  = {:10.6f}".format(2 * indent, self.hw_max / E2C) )
        if self.hw_steps != None:
            print( "{} hw_steps     = {}".format(2 * indent, int(self.hw_steps)) )
        else:   
            print( "{} resolution   = {}".format(2 * indent, int(self.resolution)) )
        print( "{} cal_method   = {}".format(2 * indent, self.cal_method))
        #condition of integral
        if self.cal_method == "integral":
            print( "{} limit (s)    = {:10.6e}".format(2 * indent, self.limit) )
            print( "{} tolerance    = {:10.6e}".format(2 * indent, self.tolerance) )
        print( "{} PL_hw_xim = {:10.6e}".format(2 * indent, self.PL_hw_xmin) )
        print( "{} PL_hw_xmax = {:10.6e}".format(2 * indent, self.PL_hw_xmax) )
        print( "{} S_hw_xmin = {:10.6e}".format(2 * indent, self.S_hw_xmin) )
        print( "{} S_hw_xmax = {:10.6e}".format(2 * indent, self.S_hw_xmax) )
        print()


def read_input():
    """
    read input file
    return list of input file contents (or defaults if no 'filein' is found)
    """

    # create default inp
    inp = inp_class()

    # handle filein
    filein = "pl.in"
    if os.path.exists(filein):
        with open(filein, "r") as f:
            lines = f.readlines()
        # update inp
        inp = parse_input(lines, inp)
    else:
        print( "{} No '{}': Default parameters are assumed".format(indent, filein) )

    inp.print_contents()

    #return inp.contents()
    inp.gen_hw_array()
    return inp

def parse_input(lines, inp):
    """
    parse input file
    return updated inp
    """
    for line in lines:
        if ("=" in line) and ("#" not in line.split()[0]):
            arg = line.split()[0]
            if arg == "path_gs":
                inp.path_gs = str(line.split()[2])
            elif arg == "path_ex":
                inp.path_ex = str(line.split()[2])
            elif arg == "phonon_interface":
                inp.phonon_interface = str(line.split()[2])
            elif arg == "file_phonon":
                inp.file_phonon = str(line.split()[2])
            elif arg == "zpl":
                inp.zpl = float(line.split()[2])
            elif arg == "skfile":
                inp.skfile = str(line.split()[2])
            elif arg == "smear":
                inp.smear = float(line.split()[2]) * Electron2Coulomb
            elif arg == "smear_end": #([eV] --> converted [J])
                inp.smear_end = float(line.split()[2]) * Electron2Coulomb
            elif arg == "limit":
                inp.limit = float(line.split()[2])
            elif arg == "gamma":
                inp.gamma = float(line.split()[2])
            elif arg == "gamma_ev": #[eV]--> convert to Hz
                inp.gamma = float(line.split()[2]) / h_eVs
            elif arg == "tolerance":
                inp.tolerance = float(line.split()[2])
            elif arg == "hw_min":
                inp.hw_min = float(line.split()[2]) * Electron2Coulomb
            elif arg == "hw_max":
                inp.hw_max = float(line.split()[2]) * Electron2Coulomb
            elif arg == "hw_steps":
                inp.hw_steps = float(line.split()[2])
            elif arg == "resolution":
                inp.resolution= int(line.split()[2])
            elif arg == "cal_method":
                inp.cal_method = str(line.split()[2])
            elif arg == "integral_method":
                inp.integral_method = str(line.split()[2])
            elif arg == "PL_hw_xmax": #[eV] #figure of PL xrange
                inp.PL_hw_xmax = float(line.split()[2])
            elif arg == "PL_hw_xmin": #[eV] #figure of PL xrange
                inp.PL_hw_xmin = float(line.split()[2])
            elif arg == "S_hw_xmin":  #[eV] #figure of S_hw xrange
                inp.S_hw_xmin = float(line.split()[2])
            elif arg == "S_hw_xmax": #[eV] #figure of S_hw xrange
                inp.S_hw_xmax = float(line.split()[2])
            else:
                print("Argument: '%s' not recognized" % arg)
    return inp
