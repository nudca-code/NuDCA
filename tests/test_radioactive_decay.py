
import sys
from pathlib import Path
sys.path.insert(0, str(Path.cwd().resolve().parent))

import os
import time
import h5py
import numpy as np
import matplotlib.pyplot as plt

from NuDCA.NuDCA.decay_database import load_decay_database
from NuDCA.NuDCA.decay_matrix import load_decay_matrix
from NuDCA.NuDCA.radioactive_decay import RadioactiveDecayBase, RadioactiveDecay
from NuDCA.NuDCA.nuclide import Nuclide, Z_to_element
from NuDCA.NuDCA.io import get_sorted_Y, filter_Y, group_and_sum

from NuDCA.kilonovae.lightcurves import KNeLightCurve


t0 = time.perf_counter()

#-----------------------------------------------------------------------------------------------------
# traj_path = '/media/yolo/CQH/Data/WinNet_v1.1.6/CCSNe/Wang/model_9.5/traj_135937'
# traj_path = 'G:/Data/WinNet_v1.1.6/CCSNe/Wang/model_9.5/traj_135937'

# def get_Y(traj_path):
#     with h5py.File(os.path.join(traj_path, 'WinNet_data.h5'), 'r') as h5file:
#         time = np.array(h5file['snapshots/time'])
#         A_snap = np.array(h5file['snapshots/A'])
#         Z_snap = np.array(h5file['snapshots/Z'])
#         Y_snap = np.array(h5file['snapshots/Y'])
        
#     return time, Y_snap, A_snap, Z_snap


# def get_heating_rate(traj_path):
#     with h5py.File(os.path.join(traj_path, 'WinNet_data.h5'), 'r') as h5file:
#         time = np.array(h5file['energy/time'])
#         heating_rate = np.array(h5file['energy/engen_tot'])
#     return time, heating_rate



# def get_nuclide_Y(nuclide, A, Z, Y):
#     A_nuclide = Nuclide(nuclide).A
#     Z_nuclide = Nuclide(nuclide).Z
#     mask = (A==A_nuclide) & (Z==Z_nuclide)
#     Y_nuclide = np.array([Y_snap[mask][0] for Y_snap in Y])
    
#     return Y_nuclide





# decay_database = load_decay_database()
# decay_matrix = load_decay_matrix()

# index = 220
# time, Y_snap, A_snap, Z_snap = get_Y(traj_path)
# Y_ini_dict = filter_Y(Z_snap, A_snap, Y_snap[index, :])

# decay_times = time[index:]
# radioactive_decay = RadioactiveDecay(Y_ini_dict, decay_database, decay_matrix)

# t, q_dot = get_heating_rate(traj_path)
# heating_rate_table = (np.array(t[index:-1]), np.array(q_dot[index:-1]))

# q_EM = radioactive_decay.heating_rates_by_type(decay_times, 'EM')
# q_LP = radioactive_decay.heating_rates_by_type(decay_times, 'LP')
# q_HP = radioactive_decay.heating_rates_by_type(decay_times, 'HP')
# q = radioactive_decay.decay_heating_rates(decay_times)
# heating_rate_decay = (np.array(decay_times), np.array(q))


# times = np.geomspace(1.e-4, 1e2, 5000) * 86400  # days
# t1, Lbol1 = KNeLightCurve(
#         lightcurve_type = 'Luminosity',       # ['Luminosity', 'Magnitude']
#         velocity_scheme= 'Default',           # Uniform distribution
#         density_scheme= 'Kasen2017',        # ['Metzger2017', 'Kasen2017', 'Wollaeger2018']
#         thermal_scheme= 'Barnes2016_Default', # ['Barnes2016_Default', 'Barnes2016_1D', 'Barnes2016_2D']
#         heating_scheme= 'Table',       # ['Korobkin2012', 'Rosswog2024', 'Table']
#         heating_rate_data = heating_rate_table,
#         mass_ejecta = 2.e-2,  # units: M_sun
#         vel_ejecta = 0.1,     # units: c
#         opacity = 0.5,                       # ['Wu2022', 'Ekanger2023', 'Tanaka_OneVar', 'Tanaka_FourVar']
#     #   luminosity_distance   # units: Mpc
#     #   wavelength            # units: cm
#     #   **kwargs   (n_shells ...)
# )(times=times)


# t2, Lbol2 = KNeLightCurve(
#         lightcurve_type = 'Luminosity',       # ['Luminosity', 'Magnitude']
#         velocity_scheme= 'Default',           # Uniform distribution
#         density_scheme= 'Kasen2017',        # ['Metzger2017', 'Kasen2017', 'Wollaeger2018']
#         thermal_scheme= 'Barnes2016_Default', # ['Barnes2016_Default', 'Barnes2016_1D', 'Barnes2016_2D']
#         heating_scheme= 'Table',       # ['Korobkin2012', 'Rosswog2024', 'Table']
#         heating_rate_data = heating_rate_decay,
#         mass_ejecta = 2.e-2,  # units: M_sun
#         vel_ejecta = 0.1,     # units: c
#         opacity = 0.5,                       # ['Wu2022', 'Ekanger2023', 'Tanaka_OneVar', 'Tanaka_FourVar']
#     #   luminosity_distance   # units: Mpc
#     #   wavelength            # units: cm
#     #   **kwargs   (n_shells ...)
# )(times=times)


# fig, ax = plt.subplots(figsize=(8, 6), dpi=150)
# ax.plot(t1/86400, Lbol1, label='WinNet')
# ax.plot(t2/86400, Lbol2, label='Chains')
# ax.set_xlim(1e-4, 1e2)
# ax.set_ylim(1e37, 1e42)
# ax.set_xscale('log')
# ax.set_yscale('log')
# ax.set_xlabel('Time [days]')
# ax.set_ylabel('Lbol [erg / s]')

# ax.legend()
# plt.tight_layout()
# plt.show()



# fig, ax = plt.subplots(figsize=(8, 6), dpi=150)

# ax.plot(decay_times/86400, q_EM, ls='--', label=r'EM ($\gamma$-decay)')
# ax.plot(decay_times/86400, q_LP, ls='--', label=r'LP ($\beta$-decay)')
# ax.plot(decay_times/86400, q_HP, ls='--', label=r'HP ($\alpha$-decay)')
# ax.plot(decay_times/86400, q, label='Total')
# ax.plot(t/86400, q_dot, label='WinNet')
# ax.set_xlim(5e-5, 1.e2)
# ax.set_ylim(1e1, 1e16)
# ax.set_xscale('log')
# ax.set_yscale('log')
# ax.set_xlabel('Time [days]')
# ax.set_ylabel('Radioactive heating rate [erg / g / s]')

# ax.legend()
# plt.show()


#-----------------------------------------------------------------------------------------------------

traj_path_skynet = 'G:/Data/RadioactiveDecay/SkyNet/DD2_M13501350_M0'
traj_path_winnet = 'G:/Data/RadioactiveDecay/WinNet/DD2_M13501350_M0'

def get_Y_SkyNet(traj_path):
    with h5py.File(os.path.join(traj_path, 'skynet_output.h5'), 'r') as h5file:
        time = np.array(h5file['Time'])
        A_snap = np.array(h5file['A'])
        Z_snap = np.array(h5file['Z'])
        Y_snap = np.array(h5file['Y'])
        
    return time, Y_snap, A_snap, Z_snap


def get_heating_rate_SkyNet(traj_path):
    with h5py.File(os.path.join(traj_path, 'skynet_output.h5'), 'r') as h5file:
        time = np.array(h5file['Time'])
        heating_rate = np.array(h5file['HeatingRate'])
    return time, heating_rate


def get_Y_WinNet(traj_path):
    with h5py.File(os.path.join(traj_path, 'WinNet_data.h5'), 'r') as h5file:
        time = np.array(h5file['snapshots/time'])
        A_snap = np.array(h5file['snapshots/A'])
        Z_snap = np.array(h5file['snapshots/Z'])
        Y_snap = np.array(h5file['snapshots/Y'])
        
    return time, Y_snap, A_snap, Z_snap


def get_heating_rate_WinNet(traj_path):
    with h5py.File(os.path.join(traj_path, 'WinNet_data.h5'), 'r') as h5file:
        time = np.array(h5file['energy/time'])
        heating_rate_tot = np.array(h5file['energy/engen_tot'])
        heating_rate_beta = np.array(h5file['energy/engen_beta'])
        heating_rate_fiss = np.array(h5file['energy/engen_fiss'])
        
        
    return time, heating_rate_tot, heating_rate_beta, heating_rate_fiss


def get_nuclide_Y(nuclide, A, Z, Y):
    A_nuclide = Nuclide(nuclide).A
    Z_nuclide = Nuclide(nuclide).Z
    mask = (A==A_nuclide) & (Z==Z_nuclide)
    Y_nuclide = np.array([Y_snap[mask][0] for Y_snap in Y])
    
    return Y_nuclide

decay_database = load_decay_database()
decay_matrix = load_decay_matrix()

index = 545
# time, Y_snap, A_snap, Z_snap = get_Y_SkyNet(traj_path_skynet)
time, Y_snap, A_snap, Z_snap = get_Y_WinNet(traj_path_winnet)

Y_ini_dict = filter_Y(Z_snap, A_snap, Y_snap[index, :])

decay_times = time[index:] - time[index]
radioactive_decay = RadioactiveDecay(Y_ini_dict, decay_database, decay_matrix)

t_skynet, q_dot_skynet = get_heating_rate_SkyNet(traj_path_skynet)

(
    t_winnet,
    q_tot_winnet,
    q_beta_winnet,
    q_fiss_winnet
) = get_heating_rate_WinNet(traj_path_winnet)


# q_EM = radioactive_decay.heating_rates_by_type(decay_times, 'EM')
# q_LP = radioactive_decay.heating_rates_by_type(decay_times, 'LP')
# q_HP = radioactive_decay.heating_rates_by_type(decay_times, 'HP')
q = radioactive_decay.decay_heating_rates(decay_times)



# heating_rate_decay = (np.array(decay_times), np.array(q))

# heating_rate_table = (np.array(t_winnet[index:-1]), np.array(q_tot_winnet[index:-1]))
# times = np.geomspace(1.e-2, 1e2, 5000) * 86400  # days
# t1, Lbol1 = KNeLightCurve(
#         lightcurve_type = 'Luminosity',       # ['Luminosity', 'Magnitude']
#         velocity_scheme= 'Default',           # Uniform distribution
#         density_scheme= 'Metzger2017',        # ['Metzger2017', 'Kasen2017', 'Wollaeger2018']
#         thermal_scheme= 'Barnes2016_Default', # ['Barnes2016_Default', 'Barnes2016_1D', 'Barnes2016_2D']
#         heating_scheme= 'Table',       # ['Korobkin2012', 'Rosswog2024', 'Table']
#         heating_rate_data = heating_rate_table,
#         mass_ejecta = 2.e-3,  # units: M_sun
#         vel_ejecta = 0.2,     # units: c
#         opacity = 10.0,                       # ['Wu2022', 'Ekanger2023', 'Tanaka_OneVar', 'Tanaka_FourVar']
#     #   luminosity_distance   # units: Mpc
#     #   wavelength            # units: cm
#     #   **kwargs   (n_shells ...)
# )(times=times)


# t2, Lbol2 = KNeLightCurve(
#         lightcurve_type = 'Luminosity',       # ['Luminosity', 'Magnitude']
#         velocity_scheme= 'Default',           # Uniform distribution
#         density_scheme= 'Metzger2017',        # ['Metzger2017', 'Kasen2017', 'Wollaeger2018']
#         thermal_scheme= 'Barnes2016_Default', # ['Barnes2016_Default', 'Barnes2016_1D', 'Barnes2016_2D']
#         heating_scheme= 'Table',       # ['Korobkin2012', 'Rosswog2024', 'Table']
#         heating_rate_data = heating_rate_decay,
#         mass_ejecta = 2.e-3,  # units: M_sun
#         vel_ejecta = 0.2,     # units: c
#         opacity = 10.0,                       # ['Wu2022', 'Ekanger2023', 'Tanaka_OneVar', 'Tanaka_FourVar']
#     #   luminosity_distance   # units: Mpc
#     #   wavelength            # units: cm
#     #   **kwargs   (n_shells ...)
# )(times=times)



# fig, ax = plt.subplots(figsize=(8, 6), dpi=150)
# ax.plot(t1/86400, Lbol1, label='SkyNet')
# ax.plot(t2/86400, Lbol2, label='Chains')
# ax.set_xlim(1e-2, 1e2)
# ax.set_ylim(1e37, 1e42)
# ax.set_xscale('log')
# ax.set_yscale('log')
# ax.set_xlabel('Time [days]')
# ax.set_ylabel('Lbol [erg / s]')

# ax.legend()
# plt.tight_layout()
# plt.show()



fig, ax = plt.subplots(figsize=(8, 6), dpi=150)

t_decay = decay_times + time[index]

# ax.plot(t_decay/86400, q_EM, ls='--', label=r'EM ($\gamma$-decay)')
# ax.plot(t_decay/86400, q_LP, ls='--', label=r'LP ($\beta$-decay)')
# ax.plot(t_decay/86400, q_HP, ls='--', label=r'HP ($\alpha$-decay)')
ax.plot(t_decay/86400, q, ls='-', color='b', lw=3, label='Radioactive decay')
# ax.plot(t_skynet/86400, q_dot_skynet, label='SkyNet')
ax.plot(t_winnet/86400, q_tot_winnet, ls='--', lw=3, color='r', label='WinNet')
# ax.plot(t_winnet/86400, q_beta_winnet, ls=':', label='WinNet (beta)')
# ax.plot(t_winnet/86400, q_fiss_winnet, ls=':', label='WinNet (fission)')
ax.set_xlim(1.e-2, 1.e2)
ax.set_ylim(1e6, 1e14)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Time [days]')
ax.set_ylabel('Radioactive heating rate [erg / g / s]')

ax.legend(ncol=2)
plt.show()









# nuclides = ['Ni-56', 'Co-56',
#             'Ge-64', 'Ga-64',
#             'Zn-60', 'Cu-60']

# # # nuclides = ['Cu-64', 'Ni-64', 'Zn-64']
# colors = ['blue', 'orange', 'green', 'red', 'purple', 'brown']

# fig, ax = plt.subplots(figsize=(8, 6), dpi=150)
# for nuclide, color in zip(nuclides, colors):
#     Y_nuc = [result.nuclide_abundance[f'{nuclide}'] for result in results]
#     Y_net = get_nuclide_Y(nuclide, A_snap, Z_snap, Y_snap)
#     ax.plot(times, Y_nuc, ls='-', color=color, label=f'{nuclide}')
#     ax.plot(t[index:], Y_net[index:], ls='--', color=color)

# ax.set_xlim(1e0, 1.e9)
# ax.set_ylim(1e-10, 1e-2)
# ax.set_xscale('log')
# ax.set_yscale('log')
# ax.set_xlabel('Time [s]')
# ax.set_ylabel('Y')
# ax.set_title('Dashed lines: WinNet\nSolid lines: Chains')

# ax.legend()
# plt.tight_layout()
# plt.show()