import numpy as np

PI2 = 2 * np.pi

# Energy
MHZ_TO_KHZ = 1000

joule_to_khz = 1.50919E+30
hartree_to_khz = 6.57968E+12
rydberg_to_khz = 3.28984E+12
ev_to_khz = 2.41799E+11
kjpermole_to_khz = 2.50607E+09
kcalpermole_to_khz = 1.04854E+10
hz_to_khz = 1.00000E-03
khz_to_khz = 1.
mhz_to_khz = 1.00000E+03
ghz_to_khz = 1.00000E+06
thz_to_khz = 1.00000E+09
phz_to_khz = 1.00000E+12
wavenumber_to_khz = 2.99792E+07


_CONVERSION_TO_RADKHZ = {
    'ghz': ghz_to_khz * PI2,
    'radghz': ghz_to_khz,
    'mhz': mhz_to_khz * PI2,
    'radmhz': mhz_to_khz,
    'khz': khz_to_khz * PI2,
    'radkhz': khz_to_khz,
    'hz': hz_to_khz * PI2,
    'radhz': hz_to_khz,
    'joule': joule_to_khz * PI2,
}

m_to_angstrom = 1.00000E+10
bohr_to_angstrom = 5.29177E-01
a_to_angstrom = 1

_CONVERSION_TO_A = {
    'bohr': bohr_to_angstrom,
    'm': m_to_angstrom,
    'nm': m_to_angstrom * 1e-9,
    'a': a_to_angstrom
}

BOHR_TO_ANGSTROM = 5.29177E-01

HARTREE_TO_MHZ = 6579680000.0
M_TO_BOHR = 18897300000.0

ELECTRON_GYRO = -17608.597050  # rad / (ms * Gauss) or rad * kHz / G
HBAR = 1.05457172  # When everything else in rad, kHz, ms, G, A

COMPLEX_DTYPE = np.complex128

BARN_TO_BOHR2 = M_TO_BOHR ** 2 * 1E-28
EFG_CONVERSION = BARN_TO_BOHR2 * HARTREE_TO_MHZ * MHZ_TO_KHZ  # units to convert EFG

PI2 = np.pi * 2
HBAR_SI = 6.62607015e-34 / PI2 # J*s
BOHR_MAGNETON = 9.274009994E-24 
NUCLEAR_MAGNETON = 5.05078366E-27 # J * T^-1

# EFG (V/A^2) to (Hartree/Bohr_radius^2)
# V to Hartree
HARTREE_TO_eV = 27.21138624598
eV_TO_HARTREE = 1/HARTREE_TO_eV
ANGSTROM_TO_BOHR = 1/BOHR_TO_ANGSTROM 
# EFG(V/A^2) * eV_TO_HARTREE / (ANGSTROM_TO_BOHR^2) 

#eV_TO_V = 1e+18 

#γ_N/g_N = μ_N/h = 7.62259 MHz/T
NUCLEAR_MAGNETON_DIV_PLANKCONST = 5.05078366E-27 / 6.62607015e-34 * 1e-6




