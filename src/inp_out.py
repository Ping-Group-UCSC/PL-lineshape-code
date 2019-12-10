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
    path_to_qe = "."
    # zero-phonon line, units eV (converted to Joules)
    zpl = None
    # file to read previously calculated wk, Sk
    skfile = None
    # smearing for S(hw), units eV (converted to Joules)
    smear = 0.006 * Electron2Coulomb
    # limit used in integral, units of s
    limit = 3.5e-13
    # gamma factor in PL, units of s^-1
    gamma = 0.28571428571e13
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
            self.path_to_qe,
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
        print("%s path_to_qe   = %s" % (2 * indent, self.path_to_qe))
        print("%s zpl (eV)     = %s" % (2 * indent, str(self.zpl / E2C)))
        print("%s skfile       = %s" % (2 * indent, self.skfile))
        print("%s smear (eV)   = %s" % (2 * indent, str(self.smear / E2C)))
        print("%s limit (s)    = %s" % (2 * indent, str(self.limit)))
        print("%s gamma (s^-1) = %s" % (2 * indent, str(self.gamma)))
        print("%s tolerance    = %s" % (2 * indent, str(self.tolerance)))
        print("%s hw_min (eV)  = %s" % (2 * indent, str(self.hw_min / E2C)))
        print("%s hw_min (eV)  = %s" % (2 * indent, str(self.hw_max / E2C)))
        print("%s hw_steps     = %s" % (2 * indent, str(self.hw_steps)))


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
        print("No '%s': Default parameters are assumed" % filein)

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
            if arg == "path_to_qe":
                inp.path_to_qe = line.split()[2]
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
