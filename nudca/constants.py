# -*- coding: utf-8 -*-
import math
import astropy.units as unit
from astropy import constants as const


PI           = math.pi
IPI          = 1. / PI
TWO_PI       = 2. * PI
FOUR_PI      = 4. * PI

# CGS Units
E_MASS      = const.m_e.cgs.value            # Electron mass   [g]
N_MASS      = const.m_n.cgs.value            # Neutron mass    [g]
P_MASS      = const.m_n.cgs.value            # Proton mass     [g]
ATOMIC_MASS = const.u.cgs.value              # Atomic mass     [g]
NA_CGS      = const.N_A.cgs.value            # Avogadro constant  [1/mol]
DAY_CGS     = unit.day.cgs.scale             # Day   [s]
YEAR_CGS    = unit.a.cgs.scale               # Year  [s]
EV_CGS      = unit.eV.cgs.scale              # Electron Volt  [erg]

KM_CGS       = unit.km.cgs.scale             #            
G_CGS        = const.G.cgs.value             #
C_CGS        = const.c.cgs.value             # Speed of light  [cm/s]
M_SUN_CGS    = const.M_sun.cgs.value         # Solar mass  [g]
H_CGS        = const.h.cgs.value             # Planck constant  [erg*s]
K_B_CGS      = const.k_B.cgs.value           # Boltzmann constant  [erg/K]
SIGMA_SB_CGS = const.sigma_sb.cgs.value      # Stefan-Boltzmann constant  [erg/cm^2/s/K^4]
PC_CGS       = const.pc.cgs.value            # Parsec  [cm]
MPC_CGS      = unit.Mpc.cgs.scale            # Mpc
JY_CGS       = unit.Jy.cgs.scale             # Jy = 1e-23
FLUX_STD     = 3631.0 * JY_CGS               #


# SI Units
KM_SI = unit.km.si.scale
M_SUN_SI = const.M_sun.si.value
