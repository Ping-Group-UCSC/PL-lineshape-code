#!/usr/bin/env python3

import numpy as np
from numpy import fft
from scipy import integrate, special
from numpy import sqrt, exp, pi, cos, sin
import matplotlib.pyplot as plt

from chem_utils import f_Element_Symbol_to_Mass
from constant import hbar_Js, Ang2m, AMU2kg, Electron2Coulomb, h_Js, h_eVs
from io_package import read_cell_and_pos_auto
from libread import read_dynmat_mold


"""
Begin methods for calculating sk from qe output
"""


def calc_qk(vecR, list_pos1, list_pos2, list_delta_r):
    """
    calc qk = sum_ai sqrt(ma) (Re,ai - Rg,ai) dr_k,ai
    units are kg^1/2 * m (SI)
    returns array of qk's
    ___INPUT___
    list_pos1, list_pos2 : atomic position of gs and ex (Rg,ai and Re,ai)
    list_delta_r: list of phonon mode vector (normalized)
    """
    list_delta_pos = [
        atom1["pos"] - atom2["pos"] for atom1, atom2 in zip(list_pos1, list_pos2)
    ]
    # fit everything into -0.5 < x < 0.5 region
    for ix, pos in enumerate(list_delta_pos):
        for i, x in enumerate(pos):
            if abs(x) > 0.5:  # Shift out of boundary; shift back
                pos[i] = x - round(x, 0)
        # convert to angstrom with vecR, then convert to meter (SI)
        list_delta_pos[ix] = np.dot(vecR, pos) * Ang2m # unit : [m]
    list_delta = [
        {"species": atom["speciesfull"], "delta": pos}
        for atom, pos in zip(list_pos1, list_delta_pos)
    ]
    #print("list delta :") #same 
        #for delta in list_delta:
        #    print(delta["delta"], delta["species"], flush=True)
        #print()
        #import sys
        #sys.exit()

    Only_delta = np.array([ delta["delta"] for delta in list_delta ]) #unit : m
    Only_mass =np.array([ sqrt(f_Element_Symbol_to_Mass(delta["species"]) * AMU2kg) for delta in list_delta ]) #unit : kg^0.5

    qk = np.array([ calc_qk_part_fast(Only_mass, Only_delta, np.array(delta_r)) for delta_r in list_delta_r ]) #unit : kg^0.5 * m
    # previous computing qk function --> too slwow due to for loop
    #    qk = np.array(
    #        [
    #            calc_qk_part(k, list_delta, list_delta_r)
    #            for k, delta_r in enumerate(list_delta_r)
    #        ]
    #    )
    return qk #unit : kg^0.5 * m

def calc_qk_part_fast(Only_mass, Only_delta, delta_r):
    """
    sub function of calc_qk above
    """
    tmp = np.sum(Only_delta*delta_r, axis=1) #np.dot(vec1, vec2)
    tmp = Only_mass * tmp #multiply mass
    
    return np.sum(tmp) #sum

def calc_qk_part(k, list_delta, list_delta_r):
    """
    sub function of calc_qk above
    """
    qk = 0
    for ix, atom in enumerate(list_delta):
        vec1 = atom["delta"]
        vec2 = np.array(list_delta_r[k])[ix]
        # convert mass from AMU to kg (SI)
        mass = f_Element_Symbol_to_Mass(atom["species"]) * AMU2kg
        qk += sqrt(mass) * np.dot(vec1, vec2)
    return qk


def Sk(wk, qk):
    """
    calc sk = wk*qk^2/(2*hbar)
    sk is unitless: [wk] = Hz, [qk] = kg^(1/2)*m, [hbar] = J*s
    returns array of sk's
    """
    sk = wk * qk ** 2 / (2 * hbar_Js)
    # sk = [w * q ** 2 / (2 * hbar_Js) for w, q in zip(wk, qk)]
    return sk


def readSk_qe(pre_gs, pre_es, dyn_file, phonon_interface='qe'):
    """
    main call for reading sk from qe output
    ____INPUT___
    pre_gs : directory of ground state calculation 
    pre_ex : directory of excited state calculation
    dyn_file : file include phonon solution, see `read_dynmat_mold` method
    phonon_interface : tool to compute phonon, see `read_dynmat_mold` method
    __return___
    returns nat, wk, qk, sk
    """
    (vecR, list_pos_f), package = read_cell_and_pos_auto(pre_gs)
    (vecR, list_pos_i), package = read_cell_and_pos_auto(pre_es)
    # print (list_pos_f)
    nat = len(list_pos_f)
    # nmodes = 3*nat
   #wk, list_delta_r = read_dynmat_mold(nat, dyn_file, interface=phonon_interface)
    #print("input of phonon file is :",dyn_file)
    wk, list_delta_r = read_dynmat_mold(nat, dyn_file, interface=phonon_interface) #unit : rad/s
    wk = np.array(wk)
    # print(wk)
    qk = calc_qk(vecR, list_pos_i, list_pos_f, list_delta_r) #unit : kg^0.5 * m
    # print(qk)
    sk = Sk(wk, qk)
    return nat, wk, qk, sk, list_delta_r


"""
End methods for calculating sk from qe output
"""


"""
Begin generating function methods
"""


def S_t(t, smear, wk, sk):
    """
    method for calculating S(t) = sum_k( S_k G_k E_k )
    where (ws = smear/hbar):
        S_k     = partial HR factors (input)
        G_k     = exp(-i*wk*t - 0.5*(ws*t)**2)
        E_k     = 0.5*erfc(1/sqrt(2)*(i*ws*t-wk/ws))
    returns tuple: real(s_t), imag(s_t)

    Below G_k is expanded into real and imaginary parts:
        G_k     = exp(-0.5*(ws*t)**2)*(cos(wk*t)-i*sin(wk*t))
    """
    # used to simplify expressions
    w_s = smear / hbar_Js

    #example of error of s_t
    #1. array result (s_t)
    # (1.9623580147755204-1.8301469875332714j)
    #2. for loop result (s_t)
    # (1.9623580147755149-1.8301469875332717j)

    ##########################
    #1. using array sum
    ##########################
    test_t = np.array([t]).flatten()
    if len(test_t) == 1:
    
        #g_k_list = exp(-0.5 * (w_s * t) ** 2 -1j*wk *t) 
        g_k_list = exp(-0.5 * (w_s * t) ** 2) * (cos(wk * t) - 1j * sin(wk * t))
        e_k_list = 0.5 * special.erfc(1 / sqrt(2) * (1j * w_s * t - wk / w_s))
        s_t = np.sum(g_k_list*e_k_list*sk)
        #print("s_t :\n",s_t)

    else :
        #when use the quad_vec, the 't' parameter has multi value!!

        leng_t = len(t); leng_wk = len(wk)
        wk_extend = np.tile(wk,(leng_t,1))
        sk_extend = np.tile(sk,(leng_t,1))
        t_extend = t.reshape(leng_t,1)
        t_extend = np.tile(t_extend,leng_wk)

        g_k_list = exp(-0.5 * (w_s * t_extend) ** 2) * (cos(wk_extend * t_extend) - 1j * sin(wk_extend * t_extend))
        e_k_list = 0.5 * special.erfc(1 / sqrt(2) * (1j * w_s * t_extend - wk_extend / w_s))
        s_t = np.sum(g_k_list*e_k_list*sk_extend, axis=1)


    ###########################
    ##2. using for loop --> too slow
    ###########################
    #s_t = 0
    #for w_k, s_k in zip(wk, sk):
    #    g_k = exp(-0.5 * (w_s * t) ** 2) * (cos(w_k * t) - 1j * sin(w_k * t))
    #    e_k = 0.5 * special.erfc(1 / sqrt(2) * (1j * w_s * t - w_k / w_s))
    #    s_t += s_k * g_k * e_k


    #print("t",t)
    #print("t_extend",t_extend[0][0], t_extend[1][0])
    #print("s_t :\n",s_t)
    #sys.exit()
    return np.real(s_t), np.imag(s_t)


def G_t(t, smear, wk, sk, hr):
    """
    method for calculating G(t) = exp(S(t)-hr)
    here we expand into real and imaginary parts:
        g_t = exp(s_t_re-hr)*(cos(s_t_im)+i*sin(s_t_im))
    """
    s_t_re, s_t_im = S_t(t, smear, wk, sk)
    g_t = exp(s_t_re - hr) * (cos(s_t_im) + 1j * sin(s_t_im))
    return np.real(g_t), np.imag(g_t)


"""
End generating function methods
"""


"""
Begin spectral function methods
"""


def A_integrand(t, hw, smear, wk, sk, hr, gamma):
    """
    computes integrand:
        exp(-gamma*|t|)*(Re(G(t))*cos(E*t/hbar) - Im(G(t))*sin(E*t/hbar))
    The unit of gamma here should be same as E/hbar, which is [rad/s]
    but out input gamma is in unit of Hz, so gamma [rad/s] = gamma[Hz]*2pi (frequency*2pi = omega)
    """
    g_t_re, g_t_im = G_t(t, smear, wk, sk, hr)
    return exp(-gamma * abs(t)) * (
        g_t_re * cos(hw * t / hbar_Js) - g_t_im * sin(hw * t / hbar_Js)
    )


def A_integral(hw, limit, smear, wk, sk, hr, gamma, tolerance, integrate_method="quad_vec"):
    """
    method for calculating PL as:
        A(ZPL - E) = 1/2pi Int[ G(t) exp(iwt - gamma*|t|) ] dt
    where the limits of integration are -inf to +inf
    here we compute only the real part (imaginary defined in seperate function is zero):
        A(ZPL - E) = 1/2pi Int[ exp(-gamma*|t|)*(Re(G(t))*cos(E*t/hbar) - Im(G(t))*sin(E*t/hbar)) ] dt
    also the limits of -inf and +inf are replaces with 'limit' where 'limit' may be replaces with ~ 3.5E-13
    """
    if integrate_method == "quad_vec":
        def f(t): return A_integrand(t, hw, smear, wk, sk, hr, gamma)
        y, err = integrate.quad_vec(f, -limit, limit)
    elif integrate_method == "romberg":
        y = integrate.romberg(
            A_integrand, -limit, limit, args=(hw, smear, wk, sk, hr, gamma), tol=tolerance, vec_func=True
        )
    elif integrate_method == 'quad':
        y = integrate.quad(
            A_integrand, -limit, limit, args=(hw, smear, wk, sk, hr, gamma), epsabs=tolerance, epsrel=tolerance
        )
    else:
        print("integrate_method only supports quad_vec and romberg")
    return y


def A_hw(hw_array, zpl, limit, smear, wk, sk, hr, gamma, tolerance, integrate_method="quad_vec"):
    """
    method for gathering A(hw)
    returns list a1e
    """
    de_array = zpl - hw_array
    #using integration
    if integrate_method == "quad_vec":
        A_hw = A_integral(de_array, limit, smear, wk, sk, hr,
                          gamma, tolerance, integrate_method)
    elif integrate_method == "romberg":
        A_hw = np.array(
            [A_integral(de, limit, smear, wk, sk, hr, gamma,
                        tolerance, integrate_method) for de in de_array]
        )
    else:
        raise ValueError("integrate_method only supports quad_vec and romberg")
        
    pl_hw = (hw_array ** 3 / hbar_Js ** 3) * A_hw
    pl_hw_norm = (pl_hw - min(pl_hw)) / (max(pl_hw) - min(pl_hw))

    return A_hw, pl_hw, pl_hw_norm


def A_hw_FFT(hw_array, S_hw_array, zpl, smear, sk, hr, gamma, test_plot=False):
    """
    method for gathering A(hw) with FFT (scipy)
    hw_array    : unit [J]
    S_hw_array  : unit [1/J]
    sk          : unit [none]
    zpl         : unit [J]
    hr          : HR factor [non]
    gamma       : unit [Hz]
    returns list 
    """
    
    test_plot = True
    
    ################################################
    # unit convert
    #
    # we only use [eV], not [J], [rad/s], [Hz] in energy unit
    #   ,and use [1/eV], not [s] in time scale
    ################################################

    zpl = zpl / Electron2Coulomb #[eV]
    hw_array = hw_array / Electron2Coulomb #[eV]
    S_hw_array =S_hw_array * Electron2Coulomb #[1/eV]

    #gamma is in Hz unit; use E=hf to convert to eV unit
    gamma = gamma * h_eVs #[1/eV]

    ################################################
    # Info about x-axis in FFT
    ################################################
    ## ## freq domain ## ##
    #1. total_freq = max(hw_array) - min(hw_array)
    #2. nFs = len(hw_array) # num of point
    #3. Fs = total_freq/nFs = (max(hw_array)-min(hw_array))/len(hw_array)  # sampling freq, unit : [eV]
    nFs = len(hw_array)
    Fs = (max(hw_array) - min(hw_array)) / nFs #unit [eV]

    ## ## time domain ## ##
    #1. total_time = 1/Fs # = len(hw_array)/(max(hw_array)-min(hw_array)) #unit : [1/eV]
    #2. nTs = len(hw_array) # num of point
    #3. Ts = total_time / nTs = 1 / (max(hw_array)-min(hw_array)) #sampling time, unit : [1/eV]
    #0. time_range = np.linespace(-1/(2*Fs), 1/(2*Fs), nFs) 
    Ts = 1 / (max(hw_array)-min(hw_array)) #1/(Fs*nFs), unit : [1/eV]
    nTs = nFs
    #time_array = np.linspace(0, nFs/Fs, nFs) - nFs/Fs/2 # unit : [1/eV] #wrong way!!!!
    time_array = np.linspace(0, 1/(max(hw_array)/len(hw_array)), len(hw_array)) - 1/(max(hw_array)/len(hw_array))/2
    print("time_array :", min(time_array), max(time_array))


    ################################################
    #1. Fourier-transform of spectral function (S_t) --> time domain
    ################################################

    #delta_hw = Fs # unit [eV]
    #S_hw_array *= Fs
    print("S_hw :",min(S_hw_array), max(S_hw_array))
 
    S_t_array = fft.ifft(S_hw_array)       # we use ifft to get the time domain
                                                # unit : [none]
                                                # because of micro freq (d(hw))
    
    S_t_array = fft.ifftshift(S_t_array) # shift the fft result, to avoid the Nyquist freq.
    
    ##test 
    #time_array = []
    #for i in range(len(S_t_array)):
    #    time_array.append((1/resolution) * (i-len(S_t_array)/2))
    #print(time_array)

    if test_plot == True:
        print("S_t :",min(S_t_array), max(S_t_array))

        #plot the S_t
        plt.plot(time_array, S_t_array.real, label="Re")
        plt.plot(time_array, S_t_array.imag, label="Im")    
        #plt.xlim(0, 3.5e-13)
        plt.title("S_t")
        plt.savefig("test_S_t.png", dpi=300)
        plt.close()
    

    ################################################
    #2. generation function (G_t) --> time domain
    #   method for calculating G(t) = exp(S(t)-hr) --> hr = S(0) : HR factor
    #   hr : unit [none] #HR-factor
    ################################################
    G_t_array = np.exp(2*np.pi* S_t_array - hr) #unit : [none]
    
    print("G_t :",min(G_t_array), max(G_t_array))
    #
    if test_plot == True:
        #plot the G_t
        plt.plot(time_array, G_t_array.real, label="Re")
        plt.plot(time_array, G_t_array.imag, label="Im")
        #plt.xlim(0, 3.5e-13)
        plt.title("G_t")
        plt.savefig("test_G_t.png", dpi=300)
        plt.close()
    

    ################################################
    #3. compute the normal PL spectrum with FFT --> freq domain
    #   using FFT (numpy)
    # gamma : unit [eV]
    ################################################
    #gamma = 0.0003  # unit : [eV]
    print("gamma :", gamma)
    #print("gamma *time_array :", gamma*min(time_array), gamma*max(time_array))

    tmp_A_hw = G_t_array * np.exp(-gamma * np.abs(time_array))
    print("tmp_A_hw_1 :", min(tmp_A_hw), "~", max(tmp_A_hw))

    if test_plot == True:
        #plot the G_t * exp(gamma)
        plt.plot(time_array, tmp_A_hw.real, label="Re")
        plt.plot(time_array, tmp_A_hw.imag, label="Im")
        #plt.xlim(0, 3.5e-13)
        plt.title("G_t - gamma")
        plt.savefig("test_G_t_gamma.png", dpi=300)
        plt.close()

    tmp_A_hw = fft.fft(tmp_A_hw) 
    print("tmp_A_hw_2 :", min(tmp_A_hw), "~", max(tmp_A_hw))
    

    ################################################
    #4. shift the A_hw peak to the zpl 
    ################################################

    if test_plot:
        #plot the A_hw before shift
        plt.plot(hw_array, tmp_A_hw)
        plt.title("A_hw (before shift)")
        plt.savefig("test_before_A_hw.png", dpi=300)
        plt.close()
    
    print("zpl_ev :", zpl)
    A_hw = tmp_A_hw.copy()
    print("A_hw :", min(A_hw), max(A_hw))

    #shit the value
    # FFT[0] --> zpl peak
    for i in range(len(tmp_A_hw)):
        A_hw[(int(zpl/Fs)-i) % len(tmp_A_hw)] = tmp_A_hw[i]
    
    
    if test_plot:
        #plot the A_hw after shift
        plt.plot(hw_array, A_hw)
        plt.title("A_hw (after shift)")
        plt.savefig("test_A_hw.png", dpi=300)
        plt.close()
    
    ################################################
    #5. compute PL spectrum with A(hw)
    ################################################
    
    #pl_hw = (hw_array ** 3 / hbar_Js ** 3) * A_hw  
    pl_hw = []
    #in case of hw_array[0] = 0eV,
    for i in range(len(A_hw)):
        pl_hw += [A_hw[i] * ((i)*(Fs))**3] #[ev^3]

    if test_plot:
        #plot the A_hw after shift
        plt.plot(hw_array, np.array(pl_hw).real, label="Re")
        plt.plot(hw_array, np.array(pl_hw).imag, label="Im")
        plt.title("pl_re_im")
        plt.legend()
        plt.savefig("test_pl_re_im.png", dpi=300)
        plt.close()
    
    pl_hw = np.abs(pl_hw) #absolute

    print("pl_hw :", min(pl_hw), max(pl_hw))
    
    #norm
    pl_hw_norm = (pl_hw - min(pl_hw)) / (max(pl_hw) - min(pl_hw))

    return A_hw, pl_hw, pl_hw_norm

"""
End generating function methods
"""

###########################################################################

"""
Begin extra methods
"""


def sort_Sk(sk):
    """
    method for sorting sk by highest sk (with index)
    """
    check = [[i, s] for i, s in enumerate(sk)]
    sorted_sk = sorted(check, key=lambda x: x[1], reverse=True)
    #    print(sorted_sk[:10])
    return sorted_sk


def S_hw(_hw, smear, wk, sk):
    """
    S(_hw) = sum_k (Sk * exp(-(_hw - hwk)^2)/(2*smear))
    smear, wk, sk, are parameters
    _hw is a variable which passed in during integration
    returns S(_hw)
    """
    return sum([s * gaussian(_hw, hbar_Js * w, smear) for w, s in zip(wk, sk)])

def gaussian(x, a, sigma):
    """
    gaussian = 1/(sigma * sqrt(2*pi)) * exp(-0.5 * ((x-a)/sigma)**2)
    """
    return 1 / (sigma * sqrt(2 * pi)) * exp(-0.5 * ((x - a) / sigma) ** 2) #unit : [1/J] (sigma)


def cal_S_hw(hw_array, smear, wk, sk, smear_end):
    """
    S(_hw) = sum_k (Sk * exp(-(_hw - hwk)^2)/(2*smear))
    smear, wk, sk, are parameters
    _hw is a variable which passed in during integration
    returns S(_hw)
    """
    #hw_array : unit [J]
    #smear : unit [J]
    #wk : unit [rad/s] 
    #sk : unit [none]
    #leng_hw = len(hw_array)
    leng_wk = len(wk)
    leng_hw = len(hw_array)

    #extend hw array
    tmp_hw_array = np.reshape(hw_array, (leng_hw, 1))
    tmp_hw_array = np.tile(tmp_hw_array, (1, leng_wk))
    #print(tmp_hw_array.shape)
    #print(tmp_hw_array[0], tmp_hw_array[-1])
    
    #extend wk array
    wk_array = np.array(wk) * hbar_Js #unit : [J]
    wk_array = np.tile(wk_array, (leng_hw, 1)) #make 2d
    #print(wk_array.shape)
    #print(wk_array[0], wk_array[-1])
    
    #extend sk array
    sk_array = np.array(sk) 
    sk_array = np.tile(sk_array, (leng_hw, 1))

    #extend smear array
    smear_lowE = smear
    if smear_end != None:
        smear_highE = smear_end
    else:
        smear_highE = smear

    #if type(smear) == list or type(smear) == np.array:
    #    if len(smear) == 2:
    #        smear_highE = smear[1]
    #        smear_lowE  = smear[0]
    #    else:
    #        smear_highE = 1.5 /1000 * Electron2Coulomb #1.5meV
    #        smear_lowE  = smear
    #else:
    #    smear_highE = 1.5/1000 * Electron2Coulomb #1.5meV
    #    smear_lowE  = smear
    

    #omega-dependent broadeing for coninuum of omega
    #smear_array = np.linspace(smear_lowE, smear_highE, leng_hw, endpoint=True)
    #smear_array = np.reshape(smear_array, (leng_hw,1))
    #smear_array = np.tile(smear_array, (1, leng_wk))

    #omega-dependent broadeing for coninuum of vibrational modes
    #smear_array = np.linspace(smear_lowE, smear_highE, leng_wk, endpoint=True)

    #smear_lowE > smear_highE
    smear_array = np.array(wk) / (max(wk) - min(wk))  * (smear_highE-smear_lowE) + smear_lowE
    smear_array = np.reshape(smear_array, (1,leng_wk))
    smear_array = np.tile(smear_array, (leng_hw, 1))

    #print("wk :", wk[0], wk[-1])
    #print("smear :", smear_lowE, smear_highE)
    


    #print(smear_array.shape)
    #print(smear_array[0], smear_array[-1])

    #compute the S_hw
    tmp_S_hw = sk_array * gaussian(tmp_hw_array, wk_array, smear_array)
    tmp_S_hw = np.sum(tmp_S_hw, axis=1)  #[1/J]
    #print(tmp_S_hw.shape)
    #print(tmp_S_hw.ndim)
    #print(sk.shape)
    
    return tmp_S_hw #unit : [1/J]
    

def gen_hw_list(hw_min, hw_max, steps):
    """
    hw_min, hw_max (eV), steps = number of points
    returns list of hw's in joules
    """
    hw_min *= Electron2Coulomb
    hw_max *= Electron2Coulomb
    dhw = (hw_max - hw_min) / steps
    return np.arange(hw_min, hw_max, dhw)


def inv_part_ratio(k, list_delta_r):
    """
    inverse partial ratio for mode k
    measures the number of atoms onto which the vibrational mode is
    localized. If, e.g., only one atom vibrates for a given mode, IPR = 1
    """
    vec = np.array(list_delta_r[k])
    vec_norm = np.linalg.norm(vec, axis=1)
    part_ratio_per_atom = vec_norm**2
    part_ratio = np.sum(part_ratio_per_atom**2)
    inv_part_ratio = 1/part_ratio
    return inv_part_ratio


"""
End extra methods
"""
