from typing import Tuple, List, Optional, Union

import numpy as np
from numba import njit
from scipy.integrate import quad
from scipy.interpolate import interp1d

from nudca.constants import C_CGS, M_SUN_CGS, FOUR_PI


class DensityProfile:
    """
    A class for modeling density profiles of ejecta in astrophysical phenomena.

    Attributes:
        density_scheme (str): The density profile scheme to use (e.g., 'Metzger2017', 'Kasen2017', 'Wollaeger2018').
        mass_ejecta (float): Total mass of the ejecta in solar masses.
        vel_min (float): Minimum velocity of the ejecta in units of the speed of light (c).
        vel_max (float): Maximum velocity of the ejecta in units of the speed of light (c).
    """
    
    def __init__(
        self,
        density_scheme: str,
        mass_ejecta: float,
        vel_min: float,
        vel_max: float
    ) -> None:
        """
        Initialize the DensityProfile object.

        Args:
            density_scheme (str): The density scheme to use.
            mass_ejecta (float): Total mass of the ejecta in solar masses.
            vel_min (float): Minimum velocity of the ejecta (fraction of the speed of light).
            vel_max (float): Maximum velocity of the ejecta (fraction of the speed of light).
        """
        self.density_scheme = density_scheme
        self.mass_ejecta = mass_ejecta
        self.vel_min = vel_min
        self.vel_max = vel_max


    def check_parameters(self):
        """
        Validate that the necessary parameters are provided.

        Raises:
            ValueError: If any required parameter is missing.
        """
        if self.mass_ejecta is None or self.vel_min is None or self.vel_max is None:
            raise ValueError('mass_ejecta, vel_min, and vel_max must be provided.')


    def density_distribution(
        self,
        times: np.ndarray,
        velocity: float
    ) -> np.ndarray:
        """
        Compute the density distribution based on the chosen density scheme.

        Args:
            times (np.ndarray): Array of times (in seconds).
            velocity (float): Velocity of the ejecta (fraction of the speed of light).

        Returns:
            np.ndarray: Density distribution at the given times and velocity.

        Raises:
            ValueError: If an unknown density scheme is provided.
        """
        self.check_parameters()

        if self.density_scheme == 'Metzger2017':
            return self._density_profile_Metzger2017(times, velocity)
        elif self.density_scheme == 'Kasen2017':
            return self._density_profile_Kasen2017(times, velocity)
        elif self.density_scheme == 'Wollaeger2018':
            return self._density_profile_Wollaeger2018(times, velocity)
        else:
            raise ValueError("Unknown density scheme")
        

    def _density_profile_Metzger2017(
        self,
        times: np.ndarray,
        velocity: float
    ) -> np.ndarray:
        """
        Compute the Metzger2017 density profile.

        Args:
            times (np.ndarray): Array of times (in seconds).
            velocity (float): Velocity of the ejecta (fraction of the speed of light).

        Returns:
            np.ndarray: Density distribution based on the Metzger2017 scheme.
        """
        mass_ejecta, R_min, R_max, radius = self._setup_constants(times, velocity)

        rho0 = (
            mass_ejecta /
            (FOUR_PI * np.power(R_min, 3) * np.log(R_max / R_min))
        )
        
        return rho0 * (radius / R_min)**-3


    def _density_profile_Kasen2017(
        self,
        times: np.ndarray,
        velocity: float
    ) -> np.ndarray:
        """
        Compute the Kasen2017 density profile.

        Args:
            times (np.ndarray): Array of times (in seconds).
            velocity (float): Velocity of the ejecta (fraction of the speed of light).

        Returns:
            np.ndarray: Density distribution based on the Kasen2017 scheme.
        """
        mass_ejecta, R_min, R_max, radius = self._setup_constants(times, velocity)

        delta, n = 1.0, 10.0
        vel_transition = (
            ((3 - n) / (3 - delta)) * 
            (
                (np.power(R_max, 3 - delta) - np.power(R_min, 3 - delta)) /
                (np.power(R_max, 3 - n) - np.power(R_min, 3 - n))
            )
        )**(1.0 / (n - delta)) / times

        zeta_rho = (
            ((3 - delta) * np.power(vel_transition * times, 3 - delta)) /
            (
                (FOUR_PI) *
                (np.power(R_max, 3 - delta) - np.power(R_min, 3 - delta))
            )
        )

        rho = np.where(
            velocity < vel_transition,
            zeta_rho * mass_ejecta / (vel_transition * times)**3 * (radius / (vel_transition * times))**-delta,
            zeta_rho * mass_ejecta / (vel_transition * times)**3 * (radius / (vel_transition * times))**-n
        )
        return rho


    def _density_profile_Wollaeger2018(
        self,
        times: np.ndarray,
        velocity: float
    ) -> np.ndarray:
        """
        Compute the Wollaeger2018 density profile.

        Args:
            times (np.ndarray): Array of times (in seconds).
            velocity (float): Velocity of the ejecta (fraction of the speed of light).

        Returns:
            np.ndarray: Density distribution based on the Wollaeger2018 scheme.
        """
        mass_ejecta, R_min, R_max, radius = self._setup_constants(velocity, times)
        
        vel_max = R_max / times
        rho0 = 315.0 * mass_ejecta * vel_max**-3 / (64.0 * np.pi) * times**-3
        
        return rho0 * (1 - (radius / R_max)**2)**3


    def _setup_constants(
        self,
        times: np.ndarray,
        velocity: float
    ) -> Tuple:
        """
        Precompute constants used in density profile calculations.

        Args:
            times (np.ndarray): Array of times (in seconds).
            velocity (float): Velocity of the ejecta (fraction of the speed of light).

        Returns:
            Tuple: Contains mass_ejecta, R_min, R_max, and radius (CGS units).
        """
        mass_ejecta = self.mass_ejecta * M_SUN_CGS
        R_min = self.vel_min * C_CGS * times
        R_max = self.vel_max * C_CGS * times
        radius = velocity * C_CGS * times
        return mass_ejecta, R_min, R_max, radius
  


class VelocityProfile:
    """
    A class for modeling velocity profiles of ejecta in astrophysical phenomena.

    Attributes:
        velocity_scheme (str): The velocity scheme to use (e.g., 'Default', 'Fryer2024').
        vel_ejecta (Optional[float]): Central ejecta velocity in units of the speed of light (c).
        vel_min (Optional[float]): Minimum velocity of the ejecta in units of the speed of light (c).
        vel_max (Optional[float]): Maximum velocity of the ejecta in units of the speed of light (c).
    """
    def __init__(
        self,
        velocity_scheme: str,
        vel_ejecta: Optional[float] = None,
        vel_min: Optional[float] = None,
        vel_max: Optional[float] = None
    ) -> None:
        """
        Initialize the VelocityProfile object.

        Args:
            velocity_scheme (str): The velocity scheme to use.
            vel_ejecta (Optional[float]): Central ejecta velocity (fraction of the speed of light).
            vel_min (Optional[float]): Minimum velocity (fraction of the speed of light).
            vel_max (Optional[float]): Maximum velocity (fraction of the speed of light).
        """
        self.velocity_scheme = velocity_scheme
        self.vel_ejecta = vel_ejecta
        
        if self.velocity_scheme == 'Default':
            self.vel_min = self.vel_ejecta - 0.05
            self.vel_max = self.vel_ejecta + 0.05
        else:
            self.vel_min = vel_min
            self.vel_max = vel_max     
        # elif self.velocity_scheme == 'Fryer2024':
        #     pass
        

    def check_parameters(self):
        """
        Validate that the necessary parameters are provided.

        Raises:
            ValueError: If vel_min or vel_max is not provided.
        """
        if self.vel_min is None or self.vel_max is None:
            raise ValueError("vel_min and vel_max must be provided.")
    

    def velocity_distribution(self, n_shells: int) -> np.ndarray:
        """
        Generate a velocity distribution for the specified number of shells.

        Args:
            n_shells (int): Number of velocity shells.

        Returns:
            np.ndnarray: Array of velocities distributed between vel_min and vel_max.
        """
        self.check_parameters()
        vel_min = self.vel_min
        vel_max = self.vel_max
        vel_shells = np.linspace(vel_min, vel_max, n_shells)
        return vel_shells
    
    
    def _velocity_profile_uniform(self, vel_min: float, vel_max: float, n_shells: int) -> np.ndarray:
        """
        Generate a uniform velocity distribution.

        Args:
            vel_min (float): Minimum velocity (fraction of the speed of light).
            vel_max (float): Maximum velocity (fraction of the speed of light).
            n_shells (int): Number of velocity shells.

        Returns:
            np.ndarray: Uniformly distributed velocities.
        """
        
        return np.linspace(vel_min, vel_max, n_shells)
    
    
    def _velocity_profile_loguniform(self, vel_min: float, vel_max: float, n_shells: int) -> np.ndarray:
        """
        Generate a logarithmic-uniform velocity distribution.

        Args:
            vel_min (float): Minimum velocity (fraction of the speed of light).
            vel_max (float): Maximum velocity (fraction of the speed of light).
            n_shells (int): Number of velocity shells.

        Returns:
            np.ndarray: Logarithmically-uniformly distributed velocities.
        """
        
        return np.geomspace(vel_min, vel_max, n_shells)
    




class Geometry:
    """
    A class for modeling the geometry of ejecta, including velocity and density profiles.

    Attributes:
        velocity_profile (VelocityProfile): Instance of the VelocityProfile class.
        density_profile (DensityProfile): Instance of the DensityProfile class.
    """
    def __init__(
        self,
        velocity_scheme: str,
        density_scheme: str,
        mass_ejecta: Optional[float] = None,
        vel_ejecta: Optional[float] = None,
        vel_min: Optional[float] = None,
        vel_max: Optional[float] = None,
    ) -> None:
        """
        Initialize the Geometry object.

        Args:
            velocity_scheme (str): Velocity scheme for the ejecta.
            density_scheme (str): Density scheme for the ejecta.
            mass_ejecta (Optional[float]): Total ejecta mass in solar masses.
            vel_ejecta (Optional[float]): Central ejecta velocity in units of the speed of light.
            vel_min (Optional[float]): Minimum ejecta velocity in units of the speed of light.
            vel_max (Optional[float]): Maximum ejecta velocity in units of the speed of light.
        """
        self.velocity_profile = VelocityProfile(
            velocity_scheme = velocity_scheme,
            vel_ejecta = vel_ejecta,
            vel_min = vel_min,
            vel_max = vel_max,
        )
        
        self.density_profile = DensityProfile(
            density_scheme = density_scheme,
            mass_ejecta = mass_ejecta,
            vel_min = self.velocity_profile.vel_min,
            vel_max = self.velocity_profile.vel_max,
        )
        
        self.check_parameters()
        
        
    def check_parameters(self):
        """
        Validate parameters for velocity and density profiles.

        Raises:
            ValueError: If required parameters are missing.
        """
        self.velocity_profile.check_parameters()
        self.density_profile.check_parameters()
        

    @property
    def mass_ejecta(self) -> float:
        """Return the total ejecta mass in units of solar mass."""
        return self.density_profile.mass_ejecta

    @property
    def vel_ejecta(self) -> float:
        """Return the average ejecta velocity in units of the speed of light."""
        return self.velocity_profile.vel_ejecta
    
    @property
    def vel_min(self) -> float:
        """Return minimum ejecta velocity in units of the speed of light."""
        return self.velocity_profile.vel_min

    @property
    def vel_max(self) -> float:
        """Return maximum ejecta velocity in units of the speed of light."""
        return self.velocity_profile.vel_max
    
    
    def velocity_shells(self, n_shells: int) -> np.ndarray:
        """
        Generate velocity shells based on the velocity profile.

        Args:
            n_shells (int): Number of velocity shells.

        Returns:
            np.ndarray: Array of velocities for the shells.
        """
        return self.velocity_profile.velocity_distribution(n_shells)
    
    
    def vel_lower_shells(self, n_shells: int) -> np.ndarray:
        """
        Generate the lower bounds of velocity shells.

        Args:
            n_shells (int): Number of velocity shells.

        Returns:
            np.ndarray: Lower velocities of the shells.
        """
        return self.velocity_shells(n_shells)[:-1]
    
    
    def vel_upper_shells(self, n_shells: int) -> np.ndarray:
        """
        Generate the upper bounds of velocity shells.

        Args:
            n_shells (int): Number of velocity shells.

        Returns:
            np.ndarray: Upper velocities of the shells.
        """
        return self.velocity_shells(n_shells)[1:]
    
    
    def vel_middle_shells(self, n_shells: int) -> np.ndarray:
        return (self.vel_lower_shells(n_shells) + self.vel_upper_shells(n_shells)) / 2.0
    
    
    def density_shells(self, time, vel_shells) -> np.ndarray:
        """
        Compute the density of each shell for a given time.

        Args:
            time (float): Time since the explosion in seconds.
            vel_shells (List[float]): Velocity shells.

        Returns:
            np.ndarray: Density values for each shell.
        """
        return np.array([self.density_profile.density_distribution(time, velocity) for velocity in vel_shells])
    
    # def mass_shells(
    #     self,
    #     times: Union[float, np.ndarray],
    #     vel_lower_shells: Union[float, np.ndarray],
    #     vel_upper_shells: Union[float, np.ndarray]
    # ) -> np.ndarray:
        
    #     def integrand(r, t):
    #         vel = r / (C_CGS * t)
    #         rho = self.density_profile.density_distribution(t, vel)
    #         return 4.0 * np.pi * r**2 * rho

    #     times = np.atleast_1d(times)
    #     vel_lower_shells = np.atleast_2d(vel_lower_shells)
    #     vel_upper_shells = np.atleast_2d(vel_upper_shells)

    #     r_lower_shells = vel_lower_shells * C_CGS * times
    #     r_upper_shells = vel_upper_shells * C_CGS * times

    #     shell_mass = np.zeros((vel_lower_shells.shape[0], len(times)))

    #     for i, (r_lower_ishell, r_upper_ishell) in enumerate(zip(r_lower_shells, r_upper_shells)):
    #         for j, (time, r_lower, r_upper) in enumerate(zip(times, r_lower_ishell, r_upper_ishell)):
    #             if r_lower < r_upper:
    #                 shell_mass[i, j] = quad(integrand, r_lower, r_upper, args=(time,))[0]

    #     return shell_mass
    
    
    # def optical_depth_shells(
    #     self,
    #     times: Union[float, np.ndarray],
    #     opacity: Union[float, np.ndarray],
    #     vel_shells: Union[float, np.ndarray]
    # ) -> np.ndarray:

    #     def integrand(r, t, kappa):
    #         v = r / (C_CGS * t)
    #         rho = self.density_profile.density_distribution(t, v)
    #         return kappa * rho

    #     times = np.atleast_1d(times)
    #     vel_shells = np.atleast_2d(vel_shells)

    #     R_max = self.vel_max * C_CGS * times
    #     r_shells = vel_shells * C_CGS * times
    #     tau = np.zeros((vel_shells.shape[0], len(times)))

    #     for i, r_ishell in enumerate(r_shells):
    #         for j, (time, r_lower, r_upper) in enumerate(zip(times, r_ishell, R_max)):
    #             if r_lower < r_upper:
    #                 tau[i, j] = quad(integrand, r_lower, r_upper, args=(time, opacity,))[0]

    #     return tau
    
    
    
    def mass_shells(
        self,
        times: Union[float, np.ndarray],
        vel_lower_shells: Union[float, np.ndarray],
        vel_upper_shells: Union[float, np.ndarray]
    ) -> np.ndarray:
        times = np.atleast_1d(times)
        vel_lower_shells = np.atleast_2d(vel_lower_shells)
        vel_upper_shells = np.atleast_2d(vel_upper_shells)
        r_lower_shells = vel_lower_shells * C_CGS * times
        r_upper_shells = vel_upper_shells * C_CGS * times
        vel_shells = (vel_lower_shells + vel_upper_shells) / 2.0
        
        density = self.density_profile.density_distribution(times, vel_shells)
        volumes = (4.0 / 3.0) * np.pi * (r_upper_shells**3 - r_lower_shells**3)
        
        return density * volumes   
    
    
    def optical_depth_shells(
        self,
        times: Union[float, np.ndarray],
        opacity: Union[float, np.ndarray],
        vel_shells: Union[float, np.ndarray]
    ) -> np.ndarray:
        times = np.atleast_1d(times)
        vel_shells = np.atleast_2d(vel_shells)
        R_max = self.vel_max * C_CGS * times
        r_shells = vel_shells * C_CGS * times
        
        optical_path = R_max - r_shells
        density = self.density_profile.density_distribution(times, vel_shells)
        
        return opacity * density * optical_path
    
    
    
    def photospheric_radius(
        self,
        times: Union[float, np.ndarray],
        opacity: float,
        tau_threshold: float = 1.0,
        n_shells: int = 100
    ) -> np.ndarray:
        
        times = np.atleast_1d(times)
        vel_middle_shells = self.vel_middle_shells(n_shells)
        tau_shells = self.optical_depth_shells(times, opacity, vel_middle_shells[:, None])
        
        radius_photon = np.zeros(len(times))
        for i, time in enumerate(times):
            tau_shells_itime = tau_shells[:, i]
            if tau_shells_itime[0] > tau_threshold and tau_shells_itime[-1] < tau_threshold:
                # Interpolate to find the velocity shell where tau = tau_threshold
                vel_ishell = interp1d(tau_shells_itime, vel_middle_shells, kind="linear", fill_value="extrapolate")(tau_threshold)
            elif tau_shells_itime[-1] > tau_threshold:
                # Optical depth > optically thin criterion for all shells, use the outermost shell
                vel_ishell = vel_middle_shells[-1]
            else:
                # Optical depth < optically thin criterion for all shells, use the innermost shell
                vel_ishell = vel_middle_shells[0]

            radius_photon[i] = vel_ishell * C_CGS * time

        return radius_photon    


