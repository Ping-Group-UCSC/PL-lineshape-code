#!/usr/bin/env python
from math import pi

# Constants
Ry2eV = 13.6056917
Ha2eV = Ry2eV * 2
Ang2cm = 1e-8
Electron2Coulomb = 1.60217662e-19 #1 charge = 1.60218e-19 coulomb potental (1eV = 1.60218e-19J)
                                    # [J/eV] : [eV] -> [J]
Ha2J = Ha2eV * Electron2Coulomb
Ang2m = 1e-10 #[m/Ang] : [Ang] -> [m]
Bohr2m = 0.52917720859e-10

Ang2Bohr = Ang2m / Bohr2m
Bohr2Ang = Bohr2m / Ang2m
AMU2kg = 1.66053904e-27 #[kg/amu.] : [amu.] -> [kg]
mass_e_kg = 9.10938215e-31
AMU2me = AMU2kg / mass_e_kg
h_eVs = 4.135667662e-15 # eV*s
h_Js = 6.626070040e-34  # J*s
hbar_eVs = h_eVs / 2 / pi   # eV*s/rad
hbar_Js = h_Js / 2 / pi     # J*s/rad
Kelvin2au = 1 / 3.1577464e5
second2au = 1 / 2.418884326505e-17

THzOrd2meV = 4.13566903434377990
THzAng2meV = 4.13566903434377990 / 2 * pi
THzToCm = 33 #transfer THz to cm-1

# added
speed_of_light = 299792458.0  # [m/s]
inv_cm_to_Hz = 100 * speed_of_light # [cm/s] : [cm^-1] -> [Hz]
conv_freq_to_omega = inv_cm_to_Hz * 2 * pi # [rad*cm/s] : [cm^-1] -> [rad * Hz]
#conv_freq_to_omega = inv_cm_to_Hz / (2*pi) # Hz / rad ??

THz_to_Hz = 1e+12 # [Hz] : [THz] -> [Hz]
conv_THz_freq_to_omega = THz_to_Hz * 2 * pi # [rad*Hz] : [THz] -> [rad * Hz]
#conv_THz_freq_to_omega = THz_to_Hz / (2*pi) # Hz / rad ??
#conv_THz_freq_to_omega = THz_to_Hz 
#conversion_factor_to_THz = 15.633302 #for vasp & phonopy data
#                                     #but it doesn't need, already applied

# for formatting
indent = "    "

