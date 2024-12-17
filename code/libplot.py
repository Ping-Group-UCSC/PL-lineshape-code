#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

from constant import hbar_Js, Electron2Coulomb
from libcalc import S_hw, S_t, G_t

def plotS_hw_new(xvals, yvals, wk, sk, labelfontsize=20, ticksfontsize=14, partialHR=False, plot_unit="meV", \
            xlim = None, ylim = None):
    """
    method for plotting s_hw;
    The new method also plot partial HR's factor
    xvals : energy, unit=J
    yvals : spectral function, unit=1/J
    """
    fig, ax1 = plt.subplots()
    partialHR = True

    ## set up plot variables
    #xmin = 0
    #xmax = 200
    #dx = int(xmax / 200)

    ## generate data
    #xvals = np.array(range(xmin, xmax, dx))
    #yvals = np.array(
    #    [S_hw(hw * Electron2Coulomb * 0.001, smear, wk, sk) for hw in xvals]
    #)
    plot_unit = "eV"

    # convert to 1/meV or normalize
    if plot_unit=="eV":
        unit_xtag = "meV";   unit_ytag = "1/eV"
        add_xfactor = 1/0.001/Electron2Coulomb;    add_yfactor = Electron2Coulomb
    elif plot_unit=="meV":
        unit_xtag = "meV";      unit_ytag = "1/meV"
        add_xfactor = 1/0.001/Electron2Coulomb;  add_yfactor = 0.001 * Electron2Coulomb
        #yvals *= 1 * (Electron2Coulomb * 0.001)
    elif plot_unit=="Norm":
        unit_xtag = "meV";        unit_ytag = "Norm"
        add_xfactor = 1/0.001/Electron2Coulomb ;        add_yfactor = 1 / max(yvals)
    
    ax1.set_ylabel("$S(\hbar\omega)$ [{}]".format(unit_ytag), fontsize=labelfontsize)
    ax1.set_xlabel("$\hbar\omega$ [{}]".format(unit_xtag),fontsize=labelfontsize)

    #plot the spectral function
    ax1.plot(xvals*add_xfactor, yvals*add_yfactor, "#0044BB")

    #plot the partial HR
    if partialHR == True:
        ax2 =ax1.twinx()
        #wk : unit [rad/s]
        #wk * hbar_Js : [J]
        ax2.scatter(wk * hbar_Js * add_xfactor, sk, s=5, color='r') # x : convert rad/s to eV 
        ax2.set_ylabel('partial HR',fontsize=13)
        ax2.tick_params(axis='y', labelsize=13,rotation=90)
        
        x_lim_max = max(wk * hbar_Js * add_xfactor)
        
        
    #setting plot
    ax1.tick_params(axis='x', labelsize=ticksfontsize)
    ax1.tick_params(axis='y', labelsize=ticksfontsize,rotation=20)
    
    #x limit
    if xlim == None or xlim == (None, None):
        x_lim_max = (float(x_lim_max) // 50 +1) * 50
        ax1.set_xlim(0, x_lim_max)
    else : 
        ax1.set_xlim(xlim)

    #y limit
    if ylim == None or ylim == (None, None):
        pass
    else:   
        ax1.set_ylim(ylim)
    

    #ax1.set_xlim(min(xvals*add_xfactor), max(xvals*add_xfactor))
    #ax1.set_xlim(0, 200)

    fig.tight_layout()
    # plt.show()
    plt.savefig('S_hw.png', dpi=300)
    plt.close()

    # messy but need to dumpy this!
    dump = np.array([xvals*add_xfactor, yvals*add_yfactor]).T #unit : [meV] vs [1/eV]
    #dump = np.array([xvals, yvals]).T
    np.savetxt('S_hw.dat', dump)
    ############################################################

def plotS_hw(smear, wk, sk):
    """
    method for plotting s_hw
    """
    # set up plot variables
    xmin = 0
    xmax = 200
    dx = int(xmax / 200)
    # generate data
    xvals = np.array(range(xmin, xmax, dx))
    yvals = np.array(
        [S_hw(hw * Electron2Coulomb * 0.001, smear, wk, sk) for hw in xvals]
    )
    # convert to 1/meV or normalize
    in_meV = True
    if in_meV:
        yvals *= 1 * (Electron2Coulomb * 0.001)
    else:
        yvals *= 1 / max(yvals)
    # plot data
    plt.plot(xvals, yvals, "#0044BB")
    plt.ylabel("$S(\hbar\omega)$ (1/meV)"), plt.xlabel("$\hbar\omega$ (meV)")
    # plt.show()
    plt.savefig('S_hw.png')
    plt.close()

    # messy but need to dumpy this!
    dump = np.array([xvals, yvals]).T
    np.savetxt('S_hw.dat', dump)
    
    return None


def plotS_t(smear, wk, sk):
    """
    method for plotting s_t
    """
    # set up plot variables
    xmin = 0
    xmax = 3.5e-13
    dx = xmax / 200
    # generate data (arange allows floats)
    xvals = np.array(np.arange(xmin, xmax, dx))
    # zip(*list) unzips list of tuples
    yvals_re, yvals_im = list(zip(*[S_t(t, smear, wk, sk) for t in xvals]))
    # plot data
    plt.plot(xvals, yvals_re, "#FF0000", label="Re($S(t)$)")
    plt.plot(xvals, yvals_im, "#222222", label="Im($S(t)$)")
    plt.ylabel("$S(t)$"), plt.xlabel("$t$ (s)"), plt.legend()
    # plt.show()
    plt.savefig('S_t.png')
    plt.close()
    return None


def plotG_t(smear, wk, sk, hr):
    """
    method for plotting G_t
    """
    # set up plot variables
    xmin = 0
    xmax = 3.5e-13
    dx = xmax / 1000
    # generate data (arange allows floats)
    xvals = np.array(np.arange(xmin, xmax, dx))
    # zip(*list) unzips list of tuples
    yvals_re, yvals_im = list(zip(*[G_t(t, smear, wk, sk, hr) for t in xvals]))
    # plot data
    plt.plot(xvals, yvals_re, "#FF0000", label="Re($G(t)$)")
    plt.plot(xvals, yvals_im, "#222222", label="Im($G(t)$)")
    plt.ylabel("$G(t)$"), plt.xlabel("$t$ (s)"), plt.legend()
    # plt.show()
    plt.savefig('G_t.png')
    plt.close()
    return None


def plotA_hw(xvals, yvals):
    """
    method for plotting g_t
    """
    # convert energies to eV
    # xvals *= 1 / Electron2Coulomb
    xvals = [x / Electron2Coulomb for x in xvals]

    plt.plot(xvals, yvals, "#7434EB")
    plt.ylabel("$A(ZPL-E)$"), plt.xlabel("$E$ (eV)")
    #plt.xlim(1.5, 1.96)
    # plt.show()
    plt.savefig('pl.png')
    plt.close()
    return None

def plotA_hw_new(xvals, yvals, xlim=None,  ylim=None, labelfontsize=20, ticksfontsize=14):
    """
    method for plotting g_t
    with more selction options for xlim and ylim
    """
    # convert energies to eV
    # xvals *= 1 / Electron2Coulomb
    xvals = [x / Electron2Coulomb for x in xvals]

    plt.plot(xvals, yvals, "#7434EB")
    #plt.ylabel("$A(ZPL-E)$",fontsize=labelfontsize), plt.xlabel("$E$ [eV]",fontsize=labelfontsize)
    plt.ylabel("$PL norm(ZPL-E)$",fontsize=labelfontsize), plt.xlabel("$E$ [eV]",fontsize=labelfontsize)

    #print(np.isnan(yvals).any())

    plt.xlim(1.3, 1.8)
    #plt.xlim(1.2, 2.4)
    #plt.xlim(min(xvals), max(xvals))
    if type(xlim) != type(None) :
        plt.xlim(xlim)
    if type(ylim) != type(None) :
        plt.ylim(ylim)
    
    plt.xticks(fontsize=ticksfontsize), plt.yticks(fontsize=ticksfontsize, rotation=20)
    # plt.show()
    plt.tight_layout()
    plt.savefig('pl.png',dpi=300)
    plt.close()
    return None

def plotSimple(xvals, yvals):
    """
    generic plotting
    """
    plt.plot(xvals, yvals), plt.show()
    return None
