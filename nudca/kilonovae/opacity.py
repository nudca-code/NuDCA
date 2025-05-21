
import numpy as np
import pandas as pd
from typing import List, Union, Optional, Callable
from pathlib import Path, PurePath
from scipy.interpolate import interp1d, RegularGridInterpolator


TANAKA_DIRPATH = PurePath(Path(__file__).resolve()).parent.joinpath('Tanaka_Opacity_Data', 'Planck_Mean_Opacity.txt')


class KNOpacity(object):
    
    def __init__(self,
                 opacity_type: str,
                 opacity_method: str):
        self.opacity_type = opacity_type
        # Current Options = Tabular,
        #                   Formula
  
        self.opacity_method = opacity_method
        # Current Options = Tanaka_OneVar,
        #                   Tanaka_FourVar,
        #                   Wu2022,
        #                   Ekanger2023,
        #                   Red,
        #                   Blue,
        #                   Purple
        
    
    def __call__(self, electron_fraction=None, time=None, density=None, temperature=None, *args, **kwargs):
        return self.return_opacity(electron_fraction, time, density, temperature)

    def return_opacity(self, electron_fraction: float, time: float, density: float, temperature: float) -> float:
        if self.opacity_type.lower() == 'tabular':
            if self.opacity_method.lower() == 'tanaka_onevar':
                opacity = self.opacity_from_tabular_Tanaka_Ye(electron_fraction)
            elif self.opacity_method.lower() == 'tanaka_fourvars':
                opacity = self.opacity_from_tabular_Tanaka_Ye_t_rho_T(electron_fraction, time, density, temperature)  
            else:
                raise ValueError('Unknown opacity method\n'
                                 'Please use parameter:\n'
                                 '"Tanaka_OneVar" for Tanaka et al. 2020\n'
                                 '"Tanaka_FourVar" for Tanaka et al. 2020')
                
        elif self.opacity_type.lower() == 'formula':
            if self.opacity_method.lower() == 'wu2022':
                opacity = self.opacity_from_formula_Wu(electron_fraction)
            elif self.opacity_method.lower() == 'ekanger2023':
                opacity = self.opacity_from_formula_Ekanger(electron_fraction)
            else:
                raise ValueError('Unknown opacity method\n'
                                 'Please use parameter:\n'
                                 '"Wu2022" for Wu et al. 2022\n'
                                 '"Ekanger2023" for Ekanger et al. 2023')
        elif self.opacity_type.lower() == 'typical':
            if self.opacity_method.lower() == 'red':
                opacity = 10.0
            elif self.opacity_method.lower() == 'blue':
                opacity = 0.5
            elif self.opacity_method.lower() == 'purple':
                opacity = 3.0
            else:
                raise ValueError('Unknown opacity method\n'
                                 'Please use parameter:\n'
                                 '"Red" for Red\n'
                                 '"Blue" for Blue\n'
                                 '"Purple" for Purple')
        else:
            raise ValueError('Unknown opacity type\n'
                             'Please use parameter:\n'
                             '"Tabular" for Tabular\n'
                             '"Formula" for Formula\n'
                             '"Typical" for Typical')
        return opacity


    @staticmethod
    def opacity_from_formula_Wu(electron_fraction: Union[float, List[float], np.ndarray]) -> Union[float, List[float]]:
        '''Wu et al. 2022
            (doi: 10.1093/mnras/stac399)
        '''
        if isinstance(electron_fraction, float):
            if electron_fraction < 0. or electron_fraction > 0.5:
                raise ValueError('The electron fraction is out of range\n'
                                 + 'Please use parameter:\n'
                                 + '0. <= electron_fraction <= 0.5')

        elif isinstance(electron_fraction, (list, np.ndarray)):
            electron_fraction = np.asarray(electron_fraction)
            if np.any(electron_fraction < 0.) or np.any(electron_fraction > 0.5):
                raise ValueError('The electron fraction is out of range\n'
                                 + 'Please use parameter:\n'
                                 + '0. <= electron_fraction <= 0.5')
        else:
            raise ValueError('Wrong input type for electron_fraction\n'
                             + 'Please use parameter:\n'
                             + 'float or list or np.ndarray')
        opacity = 9 / (1 + (4*electron_fraction)**12)
        return opacity



    @staticmethod
    def opacity_from_formula_Ekanger(electron_fraction: Union[float, List[float], np.ndarray]) -> Union[float, List[float]]:
        ''' Ekanger et al. 2023
            (doi: 10.1093/mnras/stad2348)
        ''' 
        if isinstance(electron_fraction, float):
            if electron_fraction < 0. or electron_fraction > 0.5:
                raise ValueError('The electron fraction is out of range\n'
                                 + 'Please use parameter:\n'
                                 + '0. <= electron_fraction <= 0.5')

        elif isinstance(electron_fraction, (list, np.ndarray)):
            electron_fraction = np.array(electron_fraction)
            if np.any(electron_fraction < 0.) or np.any(electron_fraction > 0.5):
                raise ValueError('The electron fraction is out of range\n'
                                 + 'Please use parameter:\n'
                                 + '0 < electron_fraction <= 0.5')
        else:
            raise ValueError('Wrong input type for electron_fraction\n'
                             + 'Please use parameter:\n'
                             + 'float or list or np.ndarray')
        opacity = 211.98 * np.exp(-12.33*electron_fraction)
        return opacity
    
    
    @staticmethod
    def opacity_from_tabular_Tanaka_Ye(electron_fraction: Union[float, List[float], np.ndarray],
                                        kind: str='linear',
                                        bounds_error:Optional[bool]=None) -> Union[float, List[float]]:
        
        '''Interpolates opacity based on electron fraction using Tanaka et al. 2020 data.
    
        Args:
            electron_fraction (Union[float, List[float], np.ndarray]): 
                    Electron fraction(s) for which to calculate opacity.
            kind (str, optional): Specifies the kind of interpolation as a string 
                    ('linear', 'nearest', 'slinear', 'quadratic', 'cubic')
            bounds_error (Optional[bool], optional): 
                    If True, when interpolated values are requested outside of the domain of the input data (electron_fraction_tanaka), a ValueError is raised.
                    If False, out of bounds values are assigned fill_value.
                    Default is None.

        Returns:
            List[float]: Interpolated opacity values.
        
        Raises:
            ValueError: If electron_fraction is out of the allowed range [0, 0.5].
            ValueError: If electron_fraction is of the wrong type.
        '''
        electron_fraction_tanaka = [0.01, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50]
        opacity_tanaka           = [33.0, 32.5, 32.2, 22.30, 5.60, 5.36, 3.30, 0.96, 0.1]
        interpolator      = interp1d(electron_fraction_tanaka, opacity_tanaka,
                                    kind=kind, bounds_error=bounds_error,
                                    fill_value='extrapolate')
        
        if isinstance(electron_fraction, (float, np.floating)):
            if electron_fraction < 0. or electron_fraction > 0.5:
                raise ValueError('The electron fraction is out of range\n'
                                 + 'Please use parameter:\n'
                                 + '0.1 <= electron_fraction <= 0.5')
            else:
                opacity = interpolator(electron_fraction)
                
        elif isinstance(electron_fraction, (List, np.ndarray)):
            electron_fraction = np.asarray(electron_fraction)
            if np.any(electron_fraction < 0.) or np.any(electron_fraction > 0.5):
                raise ValueError('The electron fraction is out of range\n'
                                 + 'Please use parameter:\n'
                                 + '0 < electron_fraction <= 0.5')
            else:
                opacity = interpolator(electron_fraction).tolist()
        else:
            raise ValueError('Wrong input type for electron_fraction\n'
                             + 'Please use parameter:\n'
                             + 'float or list or np.ndarray')
        return opacity
    
    
    
    @staticmethod
    def opacity_from_tabular_Tanaka_Ye_t_rho_T(electron_fraction: float,
                                                time: float,
                                                density: float,
                                                temperature: float,
                                                method: str='linear',
                                                bounds_error: Optional[bool]=False,
                                                fill_value: Optional[float]=None) -> List[float]:
        time_tanaka = [0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
                       11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0]
        log_density_tanaka = np.arange(-19.5, -4.9, 0.5)
        temperature_tanaka = np.arange(1000, 25600, 500)
        electron_fraction_tanaka = np.arange(0.10, 0.4, 0.05)
        
        df = pd.read_csv(TANAKA_DIRPATH, sep=r'\s+', comment='#', header=None,
                         names=['electron_fraction', 'time', 'log_density', 
                                'temperature', 'kappa_bb', 'kappa_es', 'kappa_tot'])
        opacity_tanaka = df['kappa_tot'].values.reshape(len(electron_fraction_tanaka), len(time_tanaka),
                                                        len(log_density_tanaka), len(temperature_tanaka))
        
        interpolator = RegularGridInterpolator((electron_fraction_tanaka, time_tanaka, log_density_tanaka, temperature_tanaka),
                                               opacity_tanaka, method=method, bounds_error=bounds_error, fill_value=fill_value)
        
        
        if isinstance(electron_fraction, float) and isinstance(time, float) \
            and isinstance(density, (float, int)) and isinstance(temperature, (float, int)):
            
            log_density = np.log10(density)
            if electron_fraction < 0.10 or electron_fraction > 0.40 \
                or time < 0.5 or time > 20.0 \
                or log_density < -19.5 or log_density > -4.9 \
                or temperature < 1000 or temperature > 25600:
                print('The electron_fraction, time, density, temperature is out of range.\n'
                        + 'The best parameter space is:\n'
                        + '0.10 <= electron_fraction <= 0.40\n'
                        + '0.5  <= time <= 20.0 [day]\n'
                        + '-19.5 <= log_density <= -5.0 [g cm^-3]\n'
                        + '1000 <= temperature <= 25500 [K]\n'
                        + '---------------------------------------------')
                opacity = interpolator([electron_fraction, time, log_density, temperature])
                if len(opacity) == 1:
                    opacity = float(opacity[0])
                else:
                    opacity = opacity.tolist()
            else:
                opacity = interpolator([electron_fraction, time, log_density, temperature])
                if len(opacity) == 1:
                    opacity = float(opacity[0])
                else:
                    opacity = opacity.tolist()
        else:
            raise ValueError('Wrong input type for electron_fraction, time, density, temperature\n'
                             + 'Please use parameter:\n'
                             + 'float')
        return opacity
    
    
        
      
if __name__ == '__main__':
    res = KNOpacity('Formula', 'Wu2022')([0.2, 0.3])
    print(res)