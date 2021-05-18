#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

from constant import hbar_Js, Electron2Coulomb
from libcalc import S_hw, S_t, G_t


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
    plt.xlim(min(xvals), max(xvals))
    # plt.show()
    plt.savefig('pl.png')
    plt.close()
    return None


def plotSimple(xvals, yvals):
    """
    generic plotting
    """
    plt.plot(xvals, yvals), plt.show()
    return None
