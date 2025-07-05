
from typing import List, Tuple, Union, Optional

import numpy as np
from scipy.interpolate import interp1d, RegularGridInterpolator

from nudca.constants import DAY_CGS, IPI

DEFAULT_TIMES = np.linspace(1e-3, 30, 90000) * DAY_CGS


class EffectiveHeatingRate(object):
    '''
    thermal_scheme     --> str: thermal efficiency model   [Barnes_Usual, Barnes_1D, Barnes_2D, Constant]
    heating_scheme     --> str: nuclear heating model      [Korobkin2012, Kasen2017, Tabular]
    '''
    def __init__(self, 
                thermal_scheme: str,
                heating_scheme: str,
                heating_rate_data: Optional[Tuple[np.ndarray, np.ndarray]] = None,
                mass_ejecta: Optional[float] = None,
                vel_ejecta: Optional[float] = None,
                **kwargs):
        
        self.thermal_scheme = thermal_scheme
        self.heating_scheme = heating_scheme
        self.heating_rate_data = heating_rate_data
        self.mass_ejecta = mass_ejecta
        self.vel_ejecta = vel_ejecta
        self.kwargs = kwargs
    
    def __call__(self, times):
        return self.effective_heating_rate(times)

    def effective_heating_rate(
        self,
        times: Union[List[float], np.ndarray],
    ) -> Union[List[float], np.ndarray]:
        
        thermal_efficiency = ThermalizationEfficiency(thermal_scheme=self.thermal_scheme,
                                                      mass_ejecta=self.mass_ejecta,
                                                      vel_ejecta=self.vel_ejecta,
                                                      **self.kwargs)(times)
        radioactive_heating_rate = RadioactiveHeatingRate(heating_scheme=self.heating_scheme,
                                                          heating_rate_data=self.heating_rate_data,
                                                          vel_ejecta=self.vel_ejecta,
                                                          
                                                          **self.kwargs)(times)
        return thermal_efficiency * radioactive_heating_rate
    
    
    # @staticmethod
    # def log_linear_curve(x, slope, const):
    #     return slope * x + const
    
    # def rprocess_fit(self,
    #                 times: Union[List[float], np.ndarray],
    #                 mass_ejecta: float,
    #                 vel_ejecta: float,
    #                 skynet_dirpath: Optional[str]=None,
    #                 therm_efficiency: Optional[float]=None) -> list[float]:
    #     heating_rate = HeatingRate(self.heating_rate_method)(times, skynet_dirpath)
    #     thermal_efficiency = ThermalizationEfficiency(self.thermal_method)(times, mass_ejecta, vel_ejecta, therm_efficiency)
    #     TIMES_FIT = np.geomspace(min(times), max(times), 50000)
    #     rprocess_heating = thermal_efficiency * heating_rate
    #     rprocess_heating_fit = interp1d(times, rprocess_heating, kind='cubic')(TIMES_FIT)
    #     popt, _ = curve_fit(self.log_linear_curve, np.log10(TIMES_FIT), np.log10(rprocess_heating_fit))
    #     slope, const  = popt
    #     coefficient = np.power(10, const)
    #     exponent = slope
    #     return [coefficient, exponent]  




class ThermalizationEfficiency:
    def __init__(
        self,
        thermal_scheme: str = 'Barnes2016_Default',
        mass_ejecta: Optional[float]=None,
        vel_ejecta: Optional[float]=None,
        thermal_efficiency_constant: Optional[float]=None
    ) -> None:
        self.thermal_scheme = thermal_scheme
        self.mass_ejecta = mass_ejecta
        self.vel_ejecta = vel_ejecta
        self.thermal_efficiency = thermal_efficiency_constant
        
    
    def __call__(self, time):
        if self.thermal_scheme in ['Barnes2016_Default', 'Barnes2016_1D', 'Barnes2016_2D']:
            return self.thermal_efficiency_Barnes(time)
        elif self.thermal_scheme == 'Constant':
            return self.thermal_efficiency_Constant(time)
        else:
            raise ValueError("Invalid thermal_scheme.\
                              Supported schemes: 'Barnes2016_Default', 'Barnes2016_1D', 'Barnes2016_2D', 'Constant'.")


    def thermal_efficiency_coefficients(self) -> List[float]:
        if self.thermal_scheme == 'Barnes2016_Default':
            return [0.56, 0.17, 0.74]
        elif self.thermal_scheme == 'Barnes2016_1D':
            return self._interp_1d()
        elif self.thermal_scheme == 'Barnes2016_2D':
            return self._interp_2d()
        else:
            raise ValueError("Invalid thermal_scheme for coefficients. \
                             Supported scheme: 'Barnes2016_Default', 'Barnes2016_1D', 'Barnes2016_2D'.")


    def thermal_efficiency_Barnes(
        self,
        times: Union[List[float], np.ndarray]
    ) -> Union[List[float], np.ndarray]:
        a, b, d = self.thermal_efficiency_coefficients()
        times_day = times / DAY_CGS
        tmp = 2. * b * times_day**d
        eps_th = 0.36 * (np.exp(-a * times_day) + np.log(1. + tmp) / tmp)
        return eps_th
    
    
    def thermal_efficiency_Constant(
        self,
        times: Union[List[float], np.ndarray]
    ) -> Union[List[float], np.ndarray]:
        if self.thermal_efficiency is None or self.thermal_efficiency < 0. or self.thermal_efficiency > 1.:
            raise ValueError('Invalid inputs.\
                             Please use "thermal_efficiency" as input keyword argument.\
                             Thermalization efficiency must be a number inside (0, 1).')
        return np.full(len(times), self.thermal_efficiency)
    

    def _interp_1d(self):
        mass_barnes = np.asarray([1.e-3, 5.e-3, 1.e-2, 5.e-2])
        vel_barnes = np.asarray([0.1, 0.2, 0.3])
        norm = mass_barnes / vel_barnes[:, None]**2
        norm_barnes = np.sort(norm.flatten())
        a_barnes = [
            8.16, 4.52, 3.20, 2.01,
            2.19, 1.90, 1.31, 0.81,
            0.95, 0.56, 0.55, 0.27
        ]
        b_barnes = [
            1.19, 0.62, 0.45, 0.28,
            0.31, 0.28, 0.21, 0.19,
            0.15, 0.17, 0.13, 0.10
        ]
        d_barnes = [
            1.52, 1.39, 1.39, 1.12,
            1.32, 1.21, 1.13, 0.86,
            1.13, 0.74, 0.90, 0.60
        ]
        norm_value = self.mass_ejecta / self.vel_ejecta**2
        a_1d = interp1d(norm_barnes, a_barnes, kind='linear',
                        bounds_error=False, fill_value='extrapolate')(norm_value)
        b_1d = interp1d(norm_barnes, b_barnes, kind='linear',
                        bounds_error=False, fill_value='extrapolate')(norm_value)
        d_1d = interp1d(norm_barnes, d_barnes, kind='linear',
                        bounds_error=False, fill_value='extrapolate')(norm_value)
        return [a_1d, b_1d, d_1d]


    def _interp_2d(self):
        mass_barnes = [1.e-3, 5.e-3, 1.e-2, 5.e-2]
        vel_barnes = [0.1, 0.2, 0.3]
        a_barnes = [
            [2.01, 4.52, 8.16],
            [0.81, 1.90, 3.20],
            [0.56, 1.31, 2.19],
            [0.27, 0.55, 0.95]
        ]
        b_barnes = [
            [0.28, 0.62, 1.19],
            [0.19, 0.28, 0.45],
            [0.17, 0.21, 0.31],
            [0.10, 0.13, 0.15]
        ]
        d_barnes = [
            [1.12, 1.39, 1.52],
            [0.86, 1.21, 1.39],
            [0.74, 1.13, 1.32],
            [0.60, 0.90, 1.13]
        ]
        func_a_2d = RegularGridInterpolator((mass_barnes, vel_barnes), a_barnes,
                                            bounds_error=False, fill_value=None)
        func_b_2d = RegularGridInterpolator((mass_barnes, vel_barnes), b_barnes,
                                            bounds_error=False, fill_value=None)
        func_d_2d = RegularGridInterpolator((mass_barnes, vel_barnes), d_barnes,
                                            bounds_error=False, fill_value=None)
        a_2d = func_a_2d((self.mass_ejecta, self.vel_ejecta))
        b_2d = func_b_2d((self.mass_ejecta, self.vel_ejecta))
        d_2d = func_d_2d((self.mass_ejecta, self.vel_ejecta))
        return [a_2d, b_2d, d_2d]




class RadioactiveHeatingRate:
    
    def __init__(
        self,
        heating_scheme: str,
        electron_fraction: Optional[float] = None,
        vel_ejecta: Optional[float] = None,
        heating_rate_data: Optional[Tuple[np.ndarray, np.ndarray]] = None,
        **kwargs
    ) -> None:
        self.heating_scheme = heating_scheme
        self.electron_fraction = electron_fraction
        self.vel_ejecta = vel_ejecta
        self.heating_rate_data = heating_rate_data
 

    def __call__(self, times):
        return self._heating_rate(times) 
    
    
    def _heating_rate(self, times):
        if self.heating_scheme.lower() == 'korobkin2012':
            heating_rate = self._heating_rate_Korobkin2012(times)
        elif self.heating_scheme.lower() == 'rosswog2024':
            heating_rate = self._heating_rate_Rosswog2024(times)
        elif self.heating_scheme.lower() == 'kasen2017':
            heating_rate = self._heating_rate_Kasen2017(times)
        elif self.heating_scheme.lower() == 'table':
            heating_rate = self._heating_rate_Table(times)
        else:
            raise ValueError()
        return heating_rate
    
    
    def _heating_rate_Table(
        self,
        times: Union[List[float], np.ndarray],
    ) -> Union[List[float], np.ndarray]:
        
        times_table, heating_rate_table = self.heating_rate_data
        heating_rate = np.power(10, interp1d(np.log10(times_table),
                                             np.log10(heating_rate_table),
                                             kind='cubic')(np.log10(times)))
        return heating_rate
    
    
    
    @staticmethod
    def _heating_rate_Korobkin2012(
        times: Union[List[float], np.ndarray]
    ) -> Union[List[float], np.ndarray]:
        '''
        This function is a fit calculated in Korobkin et al. 2012
        (doi: 10.1111/j.1365-2966.2012.21859.x)
        '''
        eps0 = 2e18
        sigma = 0.11
        t0, alpha = 1.3, 1.3
        heating_rate = 2.0 * eps0 * (0.5 - IPI * np.arctan((times-t0)/sigma))**alpha
        return heating_rate
    
    
    @staticmethod
    def _heating_rate_Kasen2017(
        times: Union[List[float], np.ndarray]
    ) -> Union[List[float], np.ndarray]:
        alpha = 1.3
        times_day = times / DAY_CGS
        heating_rate = 1e10 * np.power(times_day, -alpha)
        return heating_rate
    
    
    def _heating_rate_Rosswog2024(
        self,
        times: Union[List[float], np.ndarray]
    ) -> Union[List[float], np.ndarray]:
        (
            eps0,
            alpha,  t0, sigma,
            alpha1, t1, sigma1, 
            logC1, tau1,
            logC2, tau2,
            logC3, tau3
        ) = self._interp_2d()
        
        eps0 = eps0 * 1.e18
        C1 = np.exp(logC1)
        C2 = np.exp(logC2)
        C3 = np.exp(logC3)
        tau1 = tau1 * 1.e3
        tau2 = tau2 * 1.e5
        tau3 = tau3 * 1.e5
        
        heating_rate = (
            eps0 * 
            (
                (0.5 - IPI * np.arctan((times - t0) / sigma))**alpha *
                (0.5 + IPI * np.arctan((times - t1) / sigma1))**alpha1
            )
            + C1 * np.exp(-times / tau1)
            + C2 * np.exp(-times / tau2)
            + C3 * np.exp(-times / tau3)
        )
        
        return heating_rate
        
        
        
        
    def _interp_2d(self):
        Ye_table = np.asarray([
            0.05, 0.10, 0.15, 0.20, 0.25,
            0.30, 0.35, 0.40, 0.45, 0.50
        ])
        vel_table = np.asarray([0.05, 0.10, 0.20, 0.30, 0.40, 0.50])
        
        eps0_table = np.asarray([
        # vej= 0.05   0.10    0.20   0.30   0.40   0.50
            [10.000, 10.000, 10.00, 10.00, 10.00, 10.00],   # Ye = 0.05
            [10.000, 10.000, 11.00, 11.00, 11.00, 11.00],   #      0.10
            [14.000, 10.000, 11.00, 11.00, 11.00, 11.00],   #      0.15
            [14.000, 10.000, 10.00, 10.00, 11.00, 11.00],   #      0.20
            [20.000, 25.000, 40.00, 38.00, 58.00, 70.00],   #      0.25
            [6.1000, 18.000, 47.10, 47.10, 74.80, 74.80],   #      0.30
            [7.3000, 7.0000, 16.30, 23.20, 43.20, 150.0],   #      0.35
            [0.0032, 0.0032, 0.008, 0.007, 0.009, 0.015],   #      0.40
            [0.2000, 0.2000, 0.600, 1.500, 1.500, 1.500],   #      0.45
            [0.4000, 1.0000, 2.000, 3.000, 3.000, 3.000]    #      0.50
        ])
        
        alpha_table = np.asarray([
            [1.370, 1.380, 1.410, 1.410, 1.410, 1.410],
            [1.410, 1.380, 1.370, 1.370, 1.370, 1.370],
            [1.410, 1.380, 1.370, 1.370, 1.370, 1.370],
            [1.360, 1.250, 1.320, 1.320, 1.340, 1.340],
            [1.440, 1.400, 1.460, 1.660, 1.600, 1.600],
            [1.360, 1.330, 1.330, 1.330, 1.374, 1.374],
            [1.400, 1.358, 1.384, 1.384, 1.384, 1.344],
            [1.800, 1.800, 2.100, 2.100, 1.900, 1.900],
            [8.000, 8.000, 7.000, 7.000, 7.000, 7.000],
            [1.400, 1.400, 1.400, 1.600, 1.600, 1.600]
        ])
          
        t0_table = np.asarray([
            [1.800, 1.400, 1.200, 1.200, 1.200, 1.200],
            [1.400, 1.000, 0.850, 0.850, 0.850, 0.850],
            [1.000, 0.800, 0.650, 0.650, 0.610, 0.610],
            [0.850, 0.600, 0.450, 0.450, 0.450, 0.450],
            [0.650, 0.380, 0.220, 0.180, 0.120, 0.095],
            [0.540, 0.310, 0.180, 0.130, 0.095, 0.081],
            [0.385, 0.235, 0.100, 0.060, 0.035, 0.025],
            [26.00, 26.00, 0.400, 0.400, 0.120, -20.0],
            [0.200, 0.120, 0.050, 0.030, 0.025, 0.021],
            [0.160, 0.080, 0.040, 0.020, 0.018, 0.016]
        ])
        
        sigma_table = np.asarray([
            [0.080, 0.080, 0.095, 0.095, 0.095, 0.095],
            [0.100, 0.080, 0.070, 0.070, 0.070, 0.070],
            [0.070, 0.080, 0.070, 0.065, 0.070, 0.070],
            [0.040, 0.030, 0.050, 0.050, 0.050, 0.050],
            [0.050, 0.030, 0.025, 0.045, 0.050, 0.050],
            [0.110, 0.040, 0.021, 0.021, 0.017, 0.017],
            [0.100, 0.094, 0.068, 0.050, 0.030, 0.010],
            [45.00, 45.00, 45.00, 45.00, 25.00, 40.00],
            [0.200, 0.120, 0.050, 0.030, 0.025, 0.021],
            [0.030, 0.015, 0.007, 0.010, 0.009, 0.007]
        ])
        
        alpha1_table = np.asarray([
            [7.500, 7.500, 7.500, 7.500, 7.500, 7.500],
            [9.000, 9.000, 7.500, 7.500, 7.000, 7.000],
            [8.000, 8.000, 7.500, 7.500, 7.000, 7.000],
            [8.000, 8.000, 7.500, 7.500, 7.000, 7.000],
            [8.000, 8.000, 5.000, 7.500, 7.000, 6.500],
            [4.500, 3.800, 4.000, 4.000, 4.000, 4.000],
            [2.400, 3.800, 3.800, 3.210, 2.910, 3.610],
            [-1.55, -1.55, -0.75, -0.75, -2.50, -5.00],
            [-1.55, -1.55, -1.55, -1.55, -1.55, -1.55],
            [3.000, 3.000, 3.000, 3.000, 3.000, 3.000]
        ])
        
        t1_table = np.asarray([
            [0.040, 0.025, 0.014, 0.010, 0.008, 0.006],
            [0.040, 0.035, 0.020, 0.012, 0.010, 0.008],
            [0.080, 0.040, 0.020, 0.012, 0.012, 0.009],
            [0.080, 0.040, 0.030, 0.018, 0.012, 0.009],
            [0.080, 0.060, 0.065, 0.028, 0.020, 0.015],
            [0.140, 0.123, 0.089, 0.060, 0.045, 0.031],
            [0.264, 0.100, 0.070, 0.055, 0.042, 0.033],
            [1.000, 1.000, 1.000, 1.000, 0.020, 0.010],
            [1.000, 1.000, 1.000, 1.000, 1.000, 1.000],
            [0.040, 0.020, 0.010, 0.002, 0.002, 0.002]
        ])
        
        sigma1_table = np.asarray([
            [0.250, 0.120, 0.045, 0.028, 0.020, 0.015],
            [0.250, 0.060, 0.035, 0.020, 0.016, 0.012],
            [0.170, 0.090, 0.035, 0.020, 0.012, 0.009],
            [0.170, 0.070, 0.035, 0.015, 0.012, 0.009],
            [0.170, 0.070, 0.050, 0.025, 0.020, 0.020],
            [0.065, 0.067, 0.053, 0.032, 0.032, 0.024],
            [0.075, 0.044, 0.030, 0.020, 0.020, 0.014],
            [10.00, 10.00, 10.00, 10.00, 0.020, 0.010],
            [10.00, 10.00, 10.00, 10.00, 10.00, 10.00],
            [0.010, 0.005, 0.002, 1.e-4, 1.e-4, 1.e-4]  
        ])
        
        logC1_table = np.asarray([
            [27.2, 27.8, 28.2, 28.2, 28.2, 28.2],
            [28.0, 27.8, 27.8, 27.8, 27.8, 27.8],
            [27.5, 27.0, 27.8, 27.8, 27.8, 27.8],
            [28.8, 28.1, 27.8, 27.8, 27.5, 27.5],
            [28.5, 28.0, 27.5, 28.5, 29.2, 29.0],
            [25.0, 27.5, 25.8, 20.9, 29.3, 1.00],
            [28.7, 27.0, 28.0, 28.0, 27.4, 25.3],
            [28.5, 29.1, 29.5, 30.1, 30.4, 29.9],
            [20.4, 20.6, 20.8, 20.9, 20.9, 21.0],
            [29.9, 30.1, 30.1, 30.2, 30.3, 30.3]
        ])
        
        tau1_table = np.asarray([
            [4.07, 4.07, 4.07, 4.07, 4.07, 4.07],
            [4.07, 4.07, 4.07, 4.07, 4.07, 4.07],
            [4.07, 4.07, 4.07, 4.07, 4.07, 4.07],
            [4.07, 4.07, 4.07, 4.07, 4.07, 4.07],
            [4.77, 4.77, 4.77, 4.77, 4.07, 4.07],
            [4.77, 4.77, 28.2, 1.03, 0.613, 1.0],
            [3.40, 14.5, 11.4, 14.3, 13.30, 13.3],
            [2.52, 2.52, 2.52, 2.52, 2.52, 2.52],
            [1.02, 1.02, 1.02, 1.02, 1.02, 1.02],
            [0.22, 0.22, 0.22, 0.22, 0.22, 0.22]
        ])
        
        logC2_table = np.asarray([
            [21.5, 21.5, 22.1, 22.1, 22.1, 22.1],
            [22.3, 21.5, 21.5, 21.8, 21.8, 21.8],
            [22.0, 21.5, 21.5, 22.0, 21.8, 21.8],
            [23.5, 22.5, 22.1, 22.0, 22.2, 22.2],
            [22.0, 22.8, 23.0, 23.0, 23.5, 23.5],
            [10.0, 0.00, 0.00, 19.8, 22.0, 21.0],
            [26.2, 14.1, 18.8, 19.1, 23.8, 19.2],
            [25.4, 25.4, 25.8, 26.0, 26.0, 25.8],
            [18.4, 18.4, 18.6, 18.6, 18.6, 18.6],
            [27.8, 28.0, 28.2, 28.2, 28.3, 28.3]
        ])
        
        tau2_table = np.asarray([
            [4.62, 4.62, 4.62, 4.62, 4.62, 4.62],
            [4.62, 4.62, 4.62, 4.62, 4.62, 4.62],
            [4.62, 4.62, 4.62, 4.62, 4.62, 4.62],
            [4.62, 4.62, 4.62, 4.62, 4.62, 4.62],
            [5.62, 5.62, 5.62, 5.62, 4.62, 4.62],
            [5.62, 5.18, 5.18, 34.7, 8.38, 22.6],
            [0.15, 4.49, 95.0, 95.0, 0.95, 146.0],
            [0.12, 0.12, 0.12, 0.12, 0.12, 0.14],
            [0.32, 0.32, 0.32, 0.32, 0.32, 0.32],
            [0.02, 0.02, 0.02, 0.02, 0.02, 0.02]
        ])
        
        logC3_table = np.asarray([
            [19.4, 19.8, 20.1, 20.1, 20.1, 20.1],
            [20.0, 19.8, 19.8, 19.8, 19.8, 19.8],
            [19.9, 19.8, 19.8, 19.8, 19.8, 19.8],
            [5.90, 9.80, 23.5, 23.5, 23.5, 23.5],
            [27.3, 26.9, 26.6, 27.4, 25.8, 25.8],
            [27.8, 26.9, 18.9, 25.4, 24.8, 25.8],
            [22.8, 17.9, 18.9, 25.4, 24.8, 25.5],
            [20.6, 20.2, 19.8, 19.2, 19.5, 18.4],
            [12.6, 13.1, 14.1, 14.5, 14.5, 14.5],
            [24.3, 24.2, 24.0, 24.0, 24.0, 23.9]
        ])
        
        tau3_table = np.asarray([
            [18.20, 18.20, 18.20, 18.20, 18.20, 18.20],
            [18.20, 18.20, 18.20, 18.20, 18.20, 18.20],
            [18.20, 18.20, 18.20, 18.20, 18.20, 18.20],
            [18.20, 18.20, 0.620, 0.620, 0.620, 0.620],
            [0.180, 0.180, 0.180, 0.180, 0.320, 0.320],
            [0.120, 0.180, 50.80, 0.180, 0.320, 0.320],
            [2.400, 51.80, 50.80, 0.180, 0.320, 0.320],
            [3.000, 2.500, 2.400, 2.400, 2.400, 60.40],
            [200.0, 200.0, 200.0, 200.0, 200.0, 200.0],
            [8.760, 8.760, 8.760, 8.760, 8.760, 8.760]
        ])
        
        params = {
            'eps0'  : eps0_table,
            'alpha' : alpha_table,
            't0'    : t0_table,
            'sigma' : sigma_table,
            'alpha1': alpha1_table,
            't1'    : t1_table,
            'sigma1': sigma1_table,
            'logC1' : logC1_table,
            'tau1'  : tau1_table,
            'logC2' : logC2_table,
            'tau2'  : tau2_table,
            'logC3' : logC3_table,
            'tau3'  : tau3_table
        }

        interpolated_values = {}
        for key, par_table in params.items():
            interpolator = RegularGridInterpolator((Ye_table, vel_table), par_table, bounds_error=False, fill_value=None)
            interpolated_values[key] = interpolator((self.electron_fraction, self.vel_ejecta))

        eps0   = interpolated_values["eps0"]
        alpha  = interpolated_values["alpha"]
        t0     = interpolated_values["t0"]
        sigma  = interpolated_values["sigma"]
        alpha1 = interpolated_values["alpha1"]
        t1     = interpolated_values["t1"]
        sigma1 = interpolated_values["sigma1"]
        logC1  = interpolated_values["logC1"]
        tau1   = interpolated_values["tau1"]
        logC2  = interpolated_values["logC2"]
        tau2   = interpolated_values["tau2"]
        logC3  = interpolated_values["logC3"]
        tau3   = interpolated_values["tau3"]

        return eps0, alpha, t0, sigma, alpha1, t1, sigma1, logC1, tau1, logC2, tau2, logC3, tau3
        
        
    

        