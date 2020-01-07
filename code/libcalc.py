#!/usr/bin/env python3

import numpy as np
import scipy
from numpy import sqrt, exp, pi, cos, sin

from chem_utils import f_Element_Symbol_to_Mass
from constant import hbar_Js, Ang2m, AMU2kg, Electron2Coulomb
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
        list_delta_pos[ix] = np.dot(vecR, pos) * Ang2m
    list_delta = [
        {"species": atom["speciesfull"], "delta": pos}
        for atom, pos in zip(list_pos1, list_delta_pos)
    ]
    qk = np.array(
        [
            calc_qk_part(k, list_delta, list_delta_r)
            for k, delta_r in enumerate(list_delta_r)
        ]
    )
    return qk


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


def readSk_qe(pre_gs, pre_es, dyn_file):
    """
    main call for reading sk from qe output
    returns nat, wk, qk, sk
    """
    (vecR, list_pos_f), package = read_cell_and_pos_auto(pre_gs)
    (vecR, list_pos_i), package = read_cell_and_pos_auto(pre_es)
    # print (list_pos_f)
    nat = len(list_pos_f)
    # nmodes = 3*nat
    wk, list_delta_r = read_dynmat_mold(nat, dyn_file)
    wk = np.array(wk)
    # print(wk)
    qk = calc_qk(vecR, list_pos_i, list_pos_f, list_delta_r)
    # print(qk)
    sk = Sk(wk, qk)
    return nat, wk, qk, sk


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
    s_t = 0
    for w_k, s_k in zip(wk, sk):
        g_k = exp(-0.5 * (w_s * t) ** 2) * (cos(w_k * t) - 1j * sin(w_k * t))
        e_k = 0.5 * scipy.special.erfc(1 / sqrt(2) * (1j * w_s * t - w_k / w_s))
        s_t += s_k * g_k * e_k
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


def A1_integrand(t, hw, smear, wk, sk, hr, gamma):
    """
    computes integrand:
        exp(-gamma*|t|)*(Re(G(t))*cos(E*t/hbar) - Im(G(t))*sin(E*t/hbar))
    """
    g_t_re, g_t_im = G_t(t, smear, wk, sk, hr)
    return exp(-gamma * abs(t)) * (
        g_t_re * cos(hw * t / hbar_Js) - g_t_im * sin(hw * t / hbar_Js)
    )


def A1_integral(hw, limit, smear, wk, sk, hr, gamma, tolerance):
    """
    method for calculating PL as:
        A(ZPL - E) = 1/2pi Int[ G(t) exp(iwt - gamma*|t|) ] dt
    where the limits of integration are -inf to +inf
    here we compute only the real part (imaginary defined in seperate function is zero):
        A(ZPL - E) = 1/2pi Int[ exp(-gamma*|t|)*(Re(G(t))*cos(E*t/hbar) - Im(G(t))*sin(E*t/hbar)) ] dt
    also the limits of -inf and +inf are replaces with 'limit' where 'limit' may be replaces with ~ 3.5E-13    
    """
    return scipy.integrate.romberg(
        A1_integrand, -limit, limit, args=(hw, smear, wk, sk, hr, gamma), tol=tolerance
    )


def A1_hw(hw_array, zpl, limit, smear, wk, sk, hr, gamma, tolerance):
    """
    method for gathering A1(E) 
    returns list a1e
    """
    de_array = zpl - hw_array
    a1_hw = np.array(
        [A1_integral(de, limit, smear, wk, sk, hr, gamma, tolerance) for de in de_array]
    )

    pl_hw = (hw_array ** 3 / hbar_Js ** 3) * a1_hw
    pl_min, pl_max = min(pl_hw), max(pl_hw)
    pl_hw_norm = (pl_hw - min(pl_hw)) / (max(pl_hw) - min(pl_hw))

    return a1_hw, pl_hw, pl_hw_norm


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
    return 1 / (sigma * sqrt(2 * pi)) * exp(-0.5 * ((x - a) / sigma) ** 2)


def gen_hw_list(hw_min, hw_max, steps):
    """
    hw_min, hw_max (eV), steps = number of points
    returns list of hw's in joules
    """
    hw_min *= Electron2Coulomb
    hw_max *= Electron2Coulomb
    dhw = (hw_max - hw_min) / steps
    return np.arange(hw_min, hw_max, dhw)


"""
End extra methods
"""
