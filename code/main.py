#!/usr/bin/env python3


import signal  # for handling ctrl-c
import sys
import os
import numpy as np
import time
import datetime

from constant import indent, Electron2Coulomb, hbar_Js
from libread import read_ZPL
from libcalc import readSk_qe, S_t, A_hw, G_t, gen_hw_list, inv_part_ratio, A_hw_FFT, cal_S_hw, S_hw
from libplot import plotS_hw, plotS_t, plotG_t, plotA_hw, plotS_hw_new, plotA_hw_new
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
    ###############################################################################
    #1. read input file
    print("\nBeginning Calculation and reading input")

    #path_gs, path_ex, phonon_interface, file_phonon, zpl, skfile, smear, limit, gamma, tolerance, hw_array = read_input() # Add
    args = read_input()
    # print(path_to_qe, zpl, skfile, smear)
    #args.integral_method = "quad_vec" #fast
    #integral_method = "romberg" #slow but accuratly #it is removed in SciPy 1.15 due to accuracy issues
    #                                           https://docs.scipy.org/doc/scipy/release/1.12.0-notes.html
    print("{}>> Computing method is '{}'\n".format(indent,args.cal_method))
    if args.cal_method == "integral":
        print("{}>> Integral method is '{}'\n".format(indent,args.integral_method))


    ###########################################################################
    ##########################################################
    #2.1. get the information about lattice and phonon (about qe)
    #the folder of ground state and excited state
    pre_gs = args.path_gs
    pre_es = args.path_ex
    if args.zpl is None:
        if os.path.exists(os.path.join(pre_gs,"OUTCAR")) and os.path.exists(os.path.join(pre_es,"OUTCAR")):
            gs_outcar = os.path.join(pre_gs,"OUTCAR")
            es_outcar = os.path.join(pre_es,"OUTCAR")
            print("Zero-phonon Line calculated from VASP output".format(indent))
            zpl = read_ZPL(gs_outcar, es_outcar, interface='vasp') * Electron2Coulomb
        else:
            print("Zero-Phonon Line calculated from QE output".format(indent))
            zpl = read_ZPL(pre_gs + ".out", pre_es + ".out", interface='qe') * Electron2Coulomb
    else:
        print("Zero-Phonon Line read from input file".format(indent))
        zpl = args.zpl * Electron2Coulomb #zpl is in J, the argument args.zpl is in eV
    print("{} ZPL = {:10.6f}\n".format(indent, (zpl / Electron2Coulomb)))

    ## 2. calc wk and sk or read it from a file
    if args.skfile == None:
        print("Calculating Sk")
        _, wk, _, sk, list_delta_r = readSk_qe(pre_gs, pre_es, args.file_phonon, phonon_interface=args.phonon_interface)
        # np.savetxt('sk.dat', np.array([wk, sk]).T)
        np.savetxt('sk.dat', np.array(
            [wk * hbar_Js / Electron2Coulomb * 1e3, sk]).T)
    else:
        if os.path.exists(args.skfile):
            print("Reading Sk")
            # enter routine to read skfile
            print("Error: routine to read skfile not yet implemented")
            pass
        else:
            raise Exception("Error: the file " +
                            str(args.skfile) + " does not exist")
    #################################################################
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


    ##################################################
    #4. calc HR factor; plot S(hw)
    #4.1 HR's factor
    hr = sum(sk)
    print("%s Huang-Rhys factor = %10.6f" % (indent, hr))
    ####################################################################
    #Compute A(hw), 
    if args.cal_method == "integral":
        #plotS_hw(args.smear, wk, sk) if plot_flag else None
        print("use integral method for A(hw)")
        if (args.S_hw_xmin == None) or (args.S_hw_xmax == None):
            xmin = 0; xmax = 2; dx = 1/1000 #resolution = 1meV
            print("plot S(hw) using default hw from 0 to 200 meV with resolusion=1/1000 eV")
        else:
            xmin = args.S_hw_xmin; xmax = args.S_hw_xmax
            resolution = 1000
            dx = 1/resolution
            print("plot S(hw) from {} to {} with resolution=1/{}".format(xmin,xmax,resolution))
        xvals_energy = np.arange(xmin, xmax, dx) #eV
        
        # hw * Electron2Coulomb --> converted to [J]
        # args.smear    --> converted to [J]
        # wk            --> rad/s
        # sk            --> no-unit
        #S_hw_array = np.array([S_hw(hw, args.smear, wk, sk) for hw in xvals_energy]) #1/J
        S_hw_array = np.array([S_hw(hw * Electron2Coulomb, args.smear, wk, sk) for hw in xvals_energy])
        #plotS_hw(args.smear, wk, sk) if plot_flag else None
        if args.plot_flag :
            plotS_hw_new(xvals_energy*Electron2Coulomb, S_hw_array, wk, sk, plot_unit="meV", xlim=None, \
                    ylim=None, \
                    labelfontsize=30, ticksfontsize=22, partialHR=args.plot_pHR) 
        # 5. S(t -> 0)
        check_hr, _ = S_t(0, args.smear, wk, sk)
        print(indent, "S(t->0) == HR?  ", hr == check_hr)
        plotS_t(args.smear, wk, sk) if args.plot_flag else None
        #6. G(t -> inf)
        print(indent, "G(t->inf) = %10.6f + %10.6f i" %
            (G_t(args.limit, args.smear, wk, sk, hr)))
        plotG_t(args.smear, wk, sk, hr) if plot_flag else None
    elif args.cal_method == "FFT":
        print("\n\tYou choose the FFT method to compute the A_hw!")
        print("\tSo, plotting of S_t and G_t is passed..\n")
        S_hw_array = cal_S_hw(args.hw_array, args.smear, wk, sk, smear_end = args.smear_end) # [J] & [1/J]
        #plotS_hw(args.smear, wk, sk) if plot_flag else None
        print("plot S(hw), hw sampling = ",args.hw_steps,"make sure you use large enough sampling for FFT")
        if (args.S_hw_xmax == None) or (args.S_hw_xmin == None):
            hw_xlim = (None,None)
        else:
            hw_xlim = (args.S_hw_xmin*1000,args.S_hw_xmax*1000) #in meV
            print("hw_xlim=",hw_xlim)
        if args.plot_flag :
            plotS_hw_new(args.hw_array, S_hw_array, wk, sk, plot_unit="meV", xlim=hw_xlim, \
                    ylim=None, \
                    labelfontsize=30, ticksfontsize=22, partialHR=args.plot_pHR) 
        #TODO: separate the G(t) and S(t) process from FFT method
        #print("S_hw_array from previous method:",S_hw_array)
        #S_hw_array = cal_S_hw(args.hw_array, args.smear, wk, sk, smear_end = args.smear_end) # [J] & [1/J]
        #print("S_hw_array",S_hw_array)
    else:
        print("\n\tError, There is no calculation method!! ",args.cal_method)
        print("\tOnly choose one method!! (integral or FFT)\n")
        sys.exit()

    ##################################################################
    #7. main routine to compute final spectral function
    print("\nCalculating A(ZPL - E)")

    if not os.path.exists("pl.dat"):
        # calculate pl and save to pl.dat
        #hw_array = gen_hw_list(1.0, 2.1, 600)
        before_time = datetime.timedelta(seconds = time.time() - start_time)
        print(indent, "time before calculate A_hw: ", before_time)
        print(indent, "now computing A_hw ...")
        ######################################
        #compute the normal PL spectrum 
        if args.cal_method == "integral":
            _, _, pl_norm = A_hw(args.hw_array, zpl, args.limit, args.smear, \
                            wk, sk, hr, args.gamma, args.tolerance, integrate_method=args.integral_method)
        elif args.cal_method == "FFT":
            _, _, pl_norm = A_hw_FFT(args.hw_array, S_hw_array, zpl, args.smear, sk, hr, args.gamma)
        ###########################################
        after_time = datetime.timedelta(seconds = time.time() - start_time)
        print(indent, "time after calculate A_hw: ", after_time, flush=True)
        with open("pl.dat", "w") as f:
            for x, y in zip(args.hw_array / Electron2Coulomb, pl_norm):
                f.write("%10.6e%s%10.6e\n" % (x, indent, y)) #unit : [eV] vs [none]

    else:
        # read data from pl.dat
        with open("pl.dat", "r") as f:
            args.hw_array = np.array(
                [float(line.split()[0]) *
                 Electron2Coulomb for line in f.readlines()]
            ) #([eV] --> converted [J])
        with open("pl.dat", "r") as f:
            pl_norm = np.array([float(line.split()[1])
                                for line in f.readlines()])

    print("\n{}PL calculation finished".format(indent))
    print("\n{}Data saved to 'pl.dat' and plot saved to 'pl.png'".format(indent))
    ###########################################################
    #8. plot the A_hw
    if  np.isnan(pl_norm).any() == True:
        print("\n{}Warning),".format(indent))
        print("{}The some element of A_hw is nan!!".format(indent))
        print("{}You should adjust the time limit in integral!!\n".format(indent))
        sys.exit()

    #plotA_hw(args.hw_array, pl_norm)
    plotA_hw_new(args.hw_array, pl_norm, xlim=(args.PL_hw_xmin, args.PL_hw_xmax), \
            labelfontsize=30, ticksfontsize=22)

    finish_time = datetime.timedelta(seconds = time.time() - start_time)
    print("\nTotal time cost for computing : ", finish_time)
    print("\nDone!")
    return None


if __name__ == "__main__":
    main()
