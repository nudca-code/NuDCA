import sys
import numpy as np
from typing import Tuple, Optional, Union
from scipy.integrate import solve_ivp

from .geometry import Geometry
from .heating_rate import EffectiveHeatingRate
from ..constants import C_CGS, H_CGS, K_B_CGS, SIGMA_SB_CGS, MPC_CGS, FLUX_STD, TWO_PI, FOUR_PI


class KNeLightCurve:
    """
    A class to model kilonova light curves based on different schemes for velocity,
    density, thermalization, and heating.

    Attributes:
        lightcurve_type (str): Specifies the type of lightcurve ('Luminosity' or 'Magnitude').
        velocity_scheme (str): The scheme for velocity profile (e.g., 'Default').
        density_scheme (str): The scheme for density distribution (e.g., 'Metzger2017').
        thermal_scheme (str): The thermalization scheme (e.g., 'Barnes2016_Default').
        heating_scheme (str): The radioactive heating scheme (e.g., 'Korobkin2012').
        mass_ejecta (float): The total ejecta mass in solar masses.
        vel_ejecta (float): The characteristic ejecta velocity.
        opacity (float): The opacity in cm^2/g.
        luminosity_distance (float): The distance to the kilonova in Mpc.
        wavelength (float): Wavelength in cm.
        n_shells (int): Number of velocity shells.
        geometry (Geometry): A Geometry object for calculating physical quantities.
    """
    
    def __init__(
        self,
        lightcurve_type: str,
        velocity_scheme: str = 'Default',
        density_scheme: str = 'Metzger2017',
        thermal_scheme: str = 'Barnes2016_Default',
        heating_scheme: str = 'Korobkin2012',
        mass_ejecta: Optional[str] = None,
        vel_ejecta: Optional[str] = None,
        opacity: Optional[float] = None,
        luminosity_distance: float = 200,
        wavelength: float = 2e-4,
        n_shells: int = 100,
        heating_rate_data: Optional[Tuple[np.ndarray, np.ndarray]] = None,
        **kwargs
    ) -> None:
        """
        Initialize the KNLightCurve class with the specified parameters.

        Args:
            lightcurve_type (str): Type of lightcurve ('Luminosity' or 'Magnitude').
            velocity_scheme (str): Velocity scheme for the model.
            density_scheme (str): Density scheme for the model.
            thermal_scheme (str): Thermalization scheme for the model.
            heating_scheme (str): Heating scheme for the model.
            mass_ejecta (Optional[str]): Ejecta mass in solar masses.
            vel_ejecta (Optional[str]): Characteristic ejecta velocity.
            opacity (Optional[float]): Opacity in cm^2/g.
            luminosity_distance (float): Distance to the source in Mpc.
            wavelength (float): Wavelength in cm.
            n_shells (int): Number of velocity shells.
        """
        self.velocity_scheme = velocity_scheme
        self.density_scheme  = density_scheme   
        self.thermal_scheme  = thermal_scheme
        self.heating_scheme  = heating_scheme
        
        self.mass_ejecta         = mass_ejecta
        self.vel_ejecta          = vel_ejecta
        self.opacity             = opacity
        self.luminosity_distance = luminosity_distance
        self.wavelength          = wavelength
        
        self.n_shells = n_shells
        
        self.heating_rate_data = heating_rate_data
        
        self.geometry = Geometry(
            velocity_scheme=self.velocity_scheme,
            density_scheme=self.density_scheme,
            mass_ejecta=self.mass_ejecta,
            vel_ejecta=self.vel_ejecta,
            
        )
        
        match lightcurve_type:
            case('Luminosity'):
                self.kilonova_emission = self.bolometric_luminosity
            case('Magnitude'):
                self.kilonova_emission = self.apparent_magnitude
            case(_):
                sys.exit('Invalid lightcurves mode\n' 
                         + 'Please use parameter:\n'
                         + '"Luminosity" to calculate bolmertic luminosity.\n'
                         + '"Magnitude" to calculate apparent magnitude')
                    
        
    def __call__(self, times):
        return self.kilonova_emission(times)
    
    
    def bolometric_luminosity(
        self,
        times: Union[float, np.ndarray]
    ) -> Tuple:
        
        vel_lower_shells = self.geometry.vel_lower_shells(self.n_shells)
        vel_upper_shells = self.geometry.vel_upper_shells(self.n_shells)
        vel_middle_shells = self.geometry.vel_middle_shells(self.n_shells)
        Ei0 = np.full(len(vel_lower_shells), 1e40)    # Initial condition

        sol = solve_ivp(
            self._diff_func,
            (times.min(), times.max()),
            Ei0,
            t_eval=times,
            args=(vel_lower_shells[:, None], vel_upper_shells[:, None],),
            vectorized=True,
        )
        
        Lbol = np.sum(
            self._luminosity_shell(
                sol.t,
                sol.y,
                vel_middle_shells[:, None],
            ),
            axis=0,
        )

        return sol.t, Lbol


    def apparent_magnitude(
        self,
        times: Union[float, np.ndarray]
    ) -> Tuple:
        
        t, Lbol = self.bolometric_luminosity(times)
        radii_photon = self.geometry.photospheric_radius(times, self.opacity)
        
        # Calculate apparent magnitude
        frequency = C_CGS / self.wavelength                        # Frequency of wavelength
        lum_distance = self.luminosity_distance * MPC_CGS          # Luminosity distance
        tmp = Lbol / (FOUR_PI * SIGMA_SB_CGS * np.multiply(radii_photon, radii_photon))
        temperature_eff = np.power(tmp, 0.25)                      #The effective Temperature
        flux_nu = (TWO_PI * H_CGS * frequency**3 * np.power(radii_photon, 2) / 
                  (C_CGS**2 * (np.exp(H_CGS * frequency / (K_B_CGS*temperature_eff))-1) * lum_distance**2))
        mag_nu = -2.5 * np.log10(flux_nu / FLUX_STD)
        
        return t, mag_nu
    

    def _luminosity_shell(
        self,
        times: Union[float, np.ndarray],
        E_shells: Union[float, np.ndarray],
        vel_shells: Union[float, np.ndarray]
    ) -> Union[float, np.ndarray]:
        
        tau_shells = self.geometry.optical_depth_shells(times, self.opacity, vel_shells)
        tlc_shells = vel_shells * times
        td_shells = tau_shells * vel_shells * times
        return E_shells / (tlc_shells + td_shells)


    def _diff_func(
        self,
        times: Union[float, np.ndarray],
        E_shells: Union[float, np.ndarray],
        vel_lower_shells: Union[float, np.ndarray],
        vel_upper_shells: Union[float, np.ndarray]
    ) -> Union[float, np.ndarray]:

        vel_shells = 0.5 * (vel_lower_shells + vel_upper_shells)
        mass_shells = self.geometry.mass_shells(times, vel_lower_shells, vel_upper_shells)
        Lbol_shells = self._luminosity_shell(
            times, E_shells, vel_shells
        )
        effective_heating_rate = EffectiveHeatingRate(
            thermal_scheme=self.thermal_scheme,
            heating_scheme=self.heating_scheme,
            mass_ejecta=self.mass_ejecta,
            vel_ejecta=self.vel_ejecta,
            heating_rate_data = self.heating_rate_data
        )(times)

        return mass_shells * effective_heating_rate - E_shells / times - Lbol_shells
