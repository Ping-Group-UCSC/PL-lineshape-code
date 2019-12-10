#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import scipy.integrate as integrate
import scipy.special as sc

# constants
h_bar = 1.054572e-34
E_zpl = 3.1162334e-19
c = 9.61305939210e-22
gamma = 0.28571428571e13


"""
read in wk and qk
"""
# qk's calculated in other py script
qk = ascii.read("../Data/new_qk_value.txt")
# inverse wavelength in cm-1 (no 2pi) these are from ph output
k = ascii.read("../Data/w_k.txt")
qk = qk["qk"]
k = k["k"]


"""
conversions and define hwk
"""

k = k * 100  # convert cm-1 to m-1
wk = 2 * np.pi * k * 299792458  # w = 2 pi c k
qk = qk * 10 ** (-10) * (1.660539 * 10 ** (-27)) ** (0.5)  # from amu to SI
hwk = h_bar * wk


"""
Read in HR factors from paper
"""

# sk2=ascii.read('../Data/2x2x2_paper_sk.csv') # unused
sk4 = ascii.read("../Data/4x4x4_paper_sk_2.csv")
# hwk2=sk4['col1']/(6.242e+18*1e+3) # unused
# shw4=ascii.read('../Data/4x4x4_shw.csv') # unused

####################################################################################################
"""
Define several functions
    not sure what Sk is for yet
    gaussian replaces dirac delta
    S_T is S(T)
    G1 is the real part of the generator function G(T)=exp(S(T)-S(O))
    h is the imaginary part of the generator function
    aa1 and aa2 are the real and imaginary parts of the spectral function
"""


def Sk(a, b):
    return (a * b ** 2) / (2 * h_bar)


def gaussian(x, mu, sig):
    norm = 1 / (np.sqrt(2 * np.pi) * sig)
    ans = norm * np.exp(-np.power(x - mu, 2.0) / (2 * np.power(sig, 2.0)))
    return ans


def S_T(s, t, d):
    step1 = (
        1
        / 2
        * np.e ** (-((t * (c ** 2 * t + 2j * d * h_bar)) / (2 * h_bar ** 2)))
        * s
        * (1 + sc.erf((-1j * c ** 2 * t + d * h_bar) / (np.sqrt(2) * c * h_bar)))
    )
    return np.sum(step1)


def G1(s, t, d):
    st = np.real(S_T(s, t, d))
    s0 = S_T(s, 0, d)
    ft = np.imag(S_T(s, t, d))
    return np.e ** (st - s0) * np.cos(ft)


def h(s, t, d):
    st = np.real(S_T(s, t, d))
    s0 = S_T(s, 0, d)
    ft = np.imag(S_T(s, t, d))
    return np.e ** (st - s0) * np.sin(ft)


def aa1(t, s, d, hw, gamma):
    g1t = G1(s, t, d)
    ht = h(s, t, d)
    return (
        (1 / (2 * np.pi))
        * np.e ** (-gamma * np.abs(t))
        * (g1t * np.cos(hw * t / h_bar) - ht * np.sin(hw * t / h_bar))
    )


def aa2(t, s, d, hw, gamma):
    g1t = G1(s, t, d)
    ht = h(s, t, d)
    return (
        (1j / (2 * np.pi))
        * np.e ** (-gamma * np.abs(t))
        * (g1t * np.sin(hw * t / h_bar) + ht * np.cos(hw * t / h_bar))
    )


####################################################################################################

"""
The calculation
"""


def calcPL(hw, ZPL, sk, hwk, fileout):
    # input hw (list, J), ZPL (float, J), sk (list, unitless), hwk (list, J)
    import time
    from constant import Electron2Coulomb

    start = time.time()
    limit = 3.5e-13
    gdat = []
    dhw = ZPL - hw
    for _dhw in dhw:
        x = integrate.romberg(
            aa1, -limit, limit, args=(sk, hwk, _dhw, gamma), tol=1e-18
        )
        gdat.append(x)
    end = time.time()
    print("   calc took %6.2f minutes" % ((end - start) / 60))

    # calc PL and rescale
    PLCg1 = (hw ** 3 / (h_bar ** 3)) * (gdat)
    PLCg1_min = min(PLCg1)
    PLCg1_max = max(PLCg1)
    PLCg1_norm = (PLCg1 - PLCg1_min) / (PLCg1_max - PLCg1_min)

    # write to output file
    if not fileout == None:
        with open(fileout, "w") as f:
            for x, y in zip(hw / Electron2Coulomb, PLCg1_norm):
                f.write("%10.6e%s%10.6e\n" % (x, "    ", y))

    return PLCg1_norm


def main():
    print()
    import time
    from constant import Electron2Coulomb
    from libplot import plotSimple

    hw2 = np.linspace(0, max(sk4["col1"] / (6.242e18 * 1e3)) * 12.5, 300)
    hw3 = E_zpl - hw2

    s = sk4["col2"]
    d = sk4["col1"] * Electron2Coulomb * 1e-3

    # run integration
    integrate_flag = False
    if integrate_flag:
        start = time.time()
        limit = 3.5e-13
        gdat = []
        for hw in hw3:
            x = integrate.romberg(aa1, -limit, limit, args=(s, d, hw, gamma), tol=1e-18)
            gdat.append(x)
        end = time.time()
        print("   calc took %6.2f minutes" % ((end - start) / 60))
    else:
        # read data from pl.dat
        with open("anikeya.dat", "r") as f:
            hw2 = [float(line.split()[0]) * Electron2Coulomb for line in f.readlines()]
        with open("anikeya.dat", "r") as f:
            PLCg1_norm = [float(line.split()[1]) for line in f.readlines()]

    if integrate_flag:
        # calc PL and rescale
        PLCg1 = (hw2 ** 3 / (h_bar ** 3)) * (gdat)
        PLCg1_min = min(PLCg1)
        PLCg1_max = max(PLCg1)
        PLCg1_norm = (PLCg1 - PLCg1_min) / (PLCg1_max - PLCg1_min)

        # write to output file
        with open("anikeya.dat", "w") as f:
            for x, y in zip(hw2 / Electron2Coulomb, PLCg1_norm):
                f.write("%10.6e%s%10.6e\n" % (x, "    ", y))

    # plot
    plotSimple(hw2 / Electron2Coulomb, PLCg1_norm)

    # for i in hw3:
    #     s = sk4["col2"]
    #     c = 9.613059e-22
    #     d = sk4["col1"] / (6.242e18 * 1e3)
    #     hw = i
    #     gamma = 0.28571428571e13 * 1
    #     x = integrate.romberg(aa1, -3.5e-13, 3.5e-13, args=(s, d, hw, gamma), tol=1e-18)
    #     gdat.append(x)
    # end = time.time()
    # print(end - start, "seconds")


if __name__ == "__main__":
    main()


#%%


"""
export final values
"""

# outfile = open("data_lumin.txt", "w")

# outfile.write("qk\n")

# for i in Qk:
#     outfile.write("%.15f" % i)
#     outfile.write("\n")

# outfile.close()


# """
# Plot
# """
# #%%
# Sk_data = ascii.read(
#     "Sk_data.csv"
# )  # alkauskas partial HR factors which were attempted to extract (note S is not reproduced)
# Shw_data = ascii.read("Shw_data.csv")  # distribution function S(hw) from alkauskas
# PL_EXP = ascii.read("PL_EXP.csv")  # experimental data from alkauskas
# PL_TH = ascii.read("PL_theory.csv")  # repeat of hse?
# PL_PBE = ascii.read("PL_TH_PBE.csv")  # alkauskas pbe lineshape
# PL_HSE = ascii.read("PL_TH_HSE.csv")  # alkauskas hse lineshape


# PLE_min = min(PL_EXP["col2"])
# PLE_max = max(PL_EXP["col2"])
# PLE_d = PL_EXP["col2"]
# PLE_norm = (PLE_d - PLE_min) / (PLE_max - PLE_min)

# PLT_min = min(PL_TH["col2"])
# PLT_max = max(PL_TH["col2"])
# PLT_d = PL_TH["col2"]
# PLT_norm = (PLT_d - PLT_min) / (PLT_max - PLT_min)

# PLB_min = min(PL_PBE["col2"])
# PLB_max = max(PL_PBE["col2"])
# PLB_d = PL_PBE["col2"]
# PLB_norm = (PLB_d - PLB_min) / (PLB_max - PLB_min)


# PLCg1 = (hw2 ** 3 / (h_bar ** 3)) * (
#     gdat
# )  # this is where you input previous calculation
# PLCg1_min = min(PLCg1)
# PLCg1_max = max(PLCg1)
# PLCg1_norm = (PLCg1 - PLCg1_min) / (PLCg1_max - PLCg1_min)


# fig = plt.figure(figsize=(20, 10))


# plt.plot(hw2 * 6.242e18, PLCg1_norm, label="PL-Cg1")
# plt.plot(PL_TH["col1"], PLT_norm, "r-", label="PL-TH")
# # plt.plot(PL_EXP['col1'],PLE_norm,'gD-',label='PL-E')
# # plt.plot(PL_PBE['col1'],PLB_norm,'ro-',label='PL_PBE')


# fig.suptitle("Photoluminescence Lineshape", fontsize=20)
# plt.xlabel("energy(eV)", fontsize=18)
# plt.ylabel(" L($\hbar$ω)(1/eV)", fontsize=16)
# plt.legend()
# plt.xlim(1.50, 2.00)

# plt.show()

