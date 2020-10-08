#!/usr/bin/env python3


import signal  # for handling ctrl-c
import sys
import os
import numpy as np
import time

from constant import indent, Electron2Coulomb, hbar_Js
from libread import read_ZPL
from libcalc import readSk_qe, S_t, A_hw, G_t, gen_hw_list, inv_part_ratio
from libplot import plotS_hw, plotS_t, plotG_t, plotA_hw
from inp_out import read_input


def handler(signal_received, frame):
    # Handle any cleanup here
    sys.stderr.write('\nSIGINT or CTRL-C detected. Exiting gracefully\n')
    sys.exit(1)


def main():
    """
    main program
    """

    # call handler if SIGINT or ctrl-c received
    signal.signal(signal.SIGINT, handler)

    plot_flag = True

    start_time = time.time()

    print("\nBeginning Calculation and reading input")
    # read input file
    path_to_qe, zpl, skfile, smear, limit, gamma, tolerance, hw_array = read_input()
    # print(path_to_qe, zpl, skfile, smear)

    pre_gs = os.path.join(path_to_qe, "relax-gs/relax")
    pre_es = os.path.join(path_to_qe, "relax-cdftup1/relax")
    dyn_file = os.path.join(path_to_qe, "relax-gs/PH/dynmat.mold")

    if zpl is None:
        print("Zero-Phonon Line calculated from QE output".format(indent))
        zpl = read_ZPL(pre_gs + ".out", pre_es + ".out") * Electron2Coulomb
    else:
        print("Zero-Phonon Line read from input file".format(indent))
        zpl = 1.945 * Electron2Coulomb
    print("{} ZPL = {:10.6f}\n".format(indent, (zpl / Electron2Coulomb)))

    # calc wk and sk or read it from a file
    if skfile is None:
        print("Calculating Sk")
        _, wk, _, sk, list_delta_r = readSk_qe(pre_gs, pre_es, dyn_file)
        # np.savetxt('sk.dat', np.array([wk, sk]).T)
        np.savetxt('sk.dat', np.array(
            [wk * hbar_Js / Electron2Coulomb * 1e3, sk]).T)
    else:
        if os.path.exists(skfile):
            print("Reading Sk")
            # enter routine to read skfile
            print("Error: routine to read skfile not yet implemented")
            pass
        else:
            raise Exception("Error: the file " +
                            str(skfile) + " does not exist")

    # calculate IPR
    modes = np.arange(0, len(list_delta_r))
    ipr = np.zeros(len(modes))
    with open("ipr.dat", "w") as ipr_file:
        ipr_file.write("# mode    phonon energy (meV)    IPR\n")
        for k in modes:
            ipr[k] = inv_part_ratio(k, list_delta_r)
            ipr_file.write("{:03d}         {:.6f}          {:.6f}\n"
                           .format(k, wk[k] * hbar_Js / Electron2Coulomb * 1e3, ipr[k]))
    print("%s Inverse Participation Ratios saved to ipr.dat" % (indent))

    # calc HR factor; plot S(hw)
    hr = sum(sk)
    print("%s Huang-Rhys factor = %10.6f" % (indent, hr))
    plotS_hw(smear, wk, sk) if plot_flag else None

    # S(t -> 0)
    check_hr, _ = S_t(0, smear, wk, sk)
    print(indent, "S(t->0) == HR?  ", hr == check_hr)
    plotS_t(smear, wk, sk) if plot_flag else None

    # G(t -> inf)
    print(indent, "G(t->inf) = %10.6f + %10.6f i" %
          (G_t(1e-12, smear, wk, sk, hr)))
    plotG_t(smear, wk, sk, hr) if plot_flag else None

    # main routine to compute final spectral function
    print("\nCalculating A(ZPL - E)")

    if not os.path.exists("pl.dat"):
        # calculate pl and save to pl.dat
        hw_array = gen_hw_list(0, 2.1, 300)
        print(indent, "time now: ", time.time() - start_time, " sec")
        _, _, pl_norm = A_hw(hw_array, zpl, limit, smear,
                             wk, sk, hr, gamma, tolerance, integrate_method="romberg")
        print(indent, "time now: ", time.time() - start_time, " sec")
        with open("pl.dat", "w") as f:
            for x, y in zip(hw_array / Electron2Coulomb, pl_norm):
                f.write("%10.6e%s%10.6e\n" % (x, indent, y))

    else:
        # read data from pl.dat
        with open("pl.dat", "r") as f:
            hw_array = np.array(
                [float(line.split()[0]) *
                 Electron2Coulomb for line in f.readlines()]
            )
        with open("pl.dat", "r") as f:
            pl_norm = np.array([float(line.split()[1])
                                for line in f.readlines()])

    print("\n{}PL calculation finished".format(indent))

    print("\n{}Data saved to 'pl.dat' and plot saved to 'pl.png'".format(indent))
    plotA_hw(hw_array, pl_norm)

    print("\nDone!")
    return None


if __name__ == "__main__":
    main()
