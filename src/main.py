#!/usr/bin/env python


import os
import numpy as np
import time

from constant import indent, Electron2Coulomb, hbar_Js
from libread import read_ZPL
from libcalc import readSk_qe, S_t, A1_hw, G_t, gen_hw_list
from libplot import plotS_hw, plotS_t, plotG_t, plotA1_hw
from inp_out import read_input

import anikeya as ak


def main():
    """
    main program
    """

    plot_flag = True

    start_time = time.time()

    print("\nBeginning Calculation and reading input")
    # read input file
    path_to_qe, zpl, skfile, smear, limit, gamma, tolerance, hw_array = read_input()
    # print(path_to_qe, zpl, skfile, smear)

    pre_gs = os.path.join(path_to_qe, "relax-gs/relax")
    pre_es = os.path.join(path_to_qe, "relax-cdftup1/relax")
    dyn_file = os.path.join(path_to_qe, "relax-gs/PH/dynmat.mold")

    if zpl == None:
        print("%s Zero-Phonon Line calculated from QE output")
        zpl = read_ZPL(pre_gs + ".out", pre_es + ".out") * Electron2Coulomb
    else:
        print("%s Zero-Phonon Line read from input file" % indent)
        zpl = 1.945 * Electron2Coulomb
    print("%s ZPL = %10.6f" % (indent, (zpl / Electron2Coulomb)))

    print()

    # calc wk and sk or read it from a file
    if skfile == None:
        print("Calculating Sk")
        _, wk, _, sk = readSk_qe(pre_gs, pre_es, dyn_file)
    else:
        if os.path.exists(skfile):
            print("Reading Sk")
            # enter routine to read skfile
            print("Error: routine to read skfile not yet implemented")
            pass
        else:
            raise Exception("Error: the file " + str(skfile) + " does not exist")

    # calc HR factor; plot S(hw)
    hr = sum(sk)
    print("%s Huang-Rhys factor = %10.6f" % (indent, hr))
    plotS_hw(smear, wk, sk) if plot_flag else None

    # S(t -> 0)
    check_hr, _ = S_t(0, smear, wk, sk)
    print(indent, "S(t->0) == HR?  ", hr == check_hr)
    plotS_t(smear, wk, sk) if plot_flag else None

    # G(t -> inf)
    print("G(t->inf) = %10.6f + %10.6f i" % (G_t(1e-12, smear, wk, sk, hr)))
    plotG_t(smear, wk, sk, hr) if plot_flag else None

    # main routine to compute final spectral function
    print("\nCalculating A(ZPL - E)")

    calc_a1_flag = True
    calc_anikeya_flag = False
    if calc_a1_flag or not os.path.exists("pl.dat"):

        hw_array = gen_hw_list(0, 2.1, 300)
        print(indent, "time now: ", time.time() - start_time, " sec")
        _, _, pl_norm = A1_hw(hw_array, zpl, limit, smear, wk, sk, hr, gamma, tolerance)
        print(indent, "time now: ", time.time() - start_time, " sec")
        with open("pl.dat", "w") as f:
            for x, y in zip(hw_array / Electron2Coulomb, pl_norm):
                f.write("%10.6e%s%10.6e\n" % (x, indent, y))

    elif calc_anikeya_flag:
        hw_array = gen_hw_list(0, 2.1, 300)
        hwk = wk * hbar_Js
        pl_norm = ak.calcPL(hw_array, zpl, sk, hwk, "ak_pl.dat")

    else:
        # read data from pl.dat
        with open("pl.dat", "r") as f:
            hw_array = np.array(
                [float(line.split()[0]) * Electron2Coulomb for line in f.readlines()]
            )
        with open("pl.dat", "r") as f:
            pl_norm = np.array([float(line.split()[1]) for line in f.readlines()])

    # if not calc_anikeya_flag:
    #     w3a1 = [hw ** 3 * a for hw, a in zip(hw_array, a1)]

    plotA1_hw(hw_array, pl_norm)

    print()
    return None


if __name__ == "__main__":
    main()

"""

At this point everything is in place to read QE and compute Sk
to-do:
    plotS_t routine
    calc g(t)
    calc a1(t)
    integrate a1(t)
    
    
    Other:
        export methods for Sk, S(t), G(t), A(E)
        function to export Sk (sorted?); significance of high Sk modes?
        read in Sk
        combine sk and wk?? may be easier and cleaner
        smear simply global i.e. defined in constants, then imported but can be changed via a method?
            -> also gamma, and hr, zpl (pros/cons to this)

"""
