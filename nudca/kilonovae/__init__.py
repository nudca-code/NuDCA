"""
Kilonovae Module for NuDCA

This module provides tools for modeling and analyzing kilonovae phenomena,
including light curve calculations, heating rates, opacity models, and
geometric considerations.

Main Features:
-------------
- Light curve calculations for kilonovae
- Heating rate models (Korobkin2012, Kasen2017, etc.)
- Thermalization efficiency calculations
- Opacity models and calculations
- Geometric modeling of kilonova ejecta

Example:
--------
>>> from nudca.kilonovae import KNeLightCurve
>>> lightcurve = KNeLightCurve(
...     lightcurve_type='Magnitude',
...     mass_ejecta=0.01,
...     vel_ejecta=0.2
... )
>>> times = np.linspace(0.1, 10, 100)  # days
>>> t, mag = lightcurve(times)
"""

from .lightcurve import KNeLightCurve
from .heating_rate import (
    EffectiveHeatingRate,
    ThermalizationEfficiency,
    RadioactiveHeatingRate
)
from .geometry import (
    Geometry,
    DensityProfile,
    VelocityProfile
)
# from .opacity import OpacityModel

__all__ = [
    'KNeLightCurve',
    'EffectiveHeatingRate',
    'ThermalizationEfficiency',
    'RadioactiveHeatingRate',
    'Geometry',
    'DensityProfile',
    'VelocityProfile',
]
