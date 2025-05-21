# Japan-Lithuania Opacity Database for Kilonova (version 1.0)
# M. Tanaka, D. Kato, G. Gaigalas, K. Kawaguchi, "Systematic opacity calculations for kilonovae" Monthly Notices of the Royal Astronomical Society 496 (2020) 1369-1392.

##################
# Opacity Table
##################

# parameter space
- Density = 10^{-19.5} - 10^{-5}     g cm^-3
  30 grid, log spacing    (Delta log(rho) = 0.5)

- Temperature = 1000 - 25500 K
  50 grid, linear spacing (Delta T = 500 K)

- Ye = 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40

- Time = 0.5, 1, 2, 3, 4 ... 19, 20 days.


# File format 
X = Ye value
Y = time

kappa_mean/: Planck Mean Opacity

- kappa_mean_yeXXX_tYY.txt
  log density (g cm^-3), temperature (K), kappa_bb_Planck (cm^2 g^-1), kappa_e_scater (cm^2 g^-1)