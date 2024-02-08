#!/usr/bin/env python3

import os
import numpy as np

from constant import Electron2Coulomb, indent


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
    # smearing for S(hw), units eV (converted to Joules)
    smear = 0.006 * Electron2Coulomb
    # limit used in integral, units of s
    limit = 3.5e-13
    # gamma factor in PL, units of s^-1
    # gamma = 0.28571428571e13
    gamma = 0.28571428571e14
    # improves accuracy of integral, unitless
    tolerance = 1e-18
    # min energy for PL spectrum
    hw_min = 0.0 * Electron2Coulomb
    # max energy for PL spectrum
    hw_max = 2.0 * Electron2Coulomb
    # number of points to produce
    hw_steps = 300

    def gen_hw_array(self):
        dhw = (self.hw_max - self.hw_min) / self.hw_steps
        return np.arange(self.hw_min, self.hw_max, dhw)

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
        )

    def print_contents(self):
        E2C = Electron2Coulomb
        print( "{} path_gs   = {}".format(2 * indent, self.path_gs) )
        print( "{} path_ex   = {}".format(2 * indent, self.path_ex) )
        print( "{} phonon_interface   = {}".format(2 * indent, self.phonon_interface) )
        print( "{} file_phonon   = {}".format(2 * indent, self.file_phonon))
        # if 'zpl' in vars():
        #     print( "{} zpl (eV)     = {}".format(2 * indent, str(self.zpl / E2C)) )
        # else:
        #     print( "{} zpl (eV)     = {}".format(2 * indent, "from input") )
        print( "{} skfile       = {}".format(2 * indent, self.skfile) )
        print( "{} smear (eV)   = {:10.6e}".format(2 * indent, self.smear / E2C) )
        print( "{} limit (s)    = {:10.6e}".format(2 * indent, self.limit) )
        print( "{} gamma (s^-1) = {:10.6e}".format(2 * indent, self.gamma) )
        print( "{} tolerance    = {:10.6e}".format(2 * indent, self.tolerance) )
        print( "{} hw_min (eV)  = {:10.6f}".format(2 * indent, self.hw_min / E2C) )
        print( "{} hw_max (eV)  = {:10.6f}".format(2 * indent, self.hw_max / E2C) )
        print( "{} hw_steps     = {}".format(2 * indent, int(self.hw_steps)) )
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

    return inp.contents()


def parse_input(lines, inp):
    """
    parse input file
    return updated inp
    """
    for line in lines:
        if ("=" in line) and ("#" not in line.split()[0]):
            arg = line.split()[0]
            if arg == "path_gs":
                inp.path_gs = line.split()[2]
            elif arg == "path_ex":
                inp.path_ex = line.split()[2]
            elif arg == "phonon_interface":
                inp.phonon_interface = line.split()[2]
            elif arg == "file_phonon":
                inp.file_phonon = line.split()[2]
            elif arg == "zpl":
                inp.zpl = float(line.split()[2]) * Electron2Coulomb
            elif arg == "skfile":
                inp.skfile = line.split()[2]
            elif arg == "smear":
                inp.smear = float(line.split()[2]) * Electron2Coulomb
            elif arg == "limit":
                inp.limit = float(line.split()[2])
            elif arg == "gamma":
                inp.gamma = float(line.split()[2])
            elif arg == "tolerance":
                inp.tolerance = float(line.split()[2])
            elif arg == "hw_min":
                inp.hw_min = float(line.split()[2]) * Electron2Coulomb
            elif arg == "hw_max":
                inp.hw_max = float(line.split()[2]) * Electron2Coulomb
            elif arg == "hw_steps":
                inp.hw_steps = float(line.split()[2])
            else:
                print("Argument: '%s' not recognized" % arg)
    return inp
