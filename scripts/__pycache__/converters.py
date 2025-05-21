import ast
import math
from abc import ABC
from typing import Union

import numpy as np
import pandas as pd
from sympy import S, Integer, nsimplify
from sympy.core.expr import Expr

from KNdecay.constants import YEAR_CGS



UNITS = {
    "ys": 1.0E-24,
    "zs": 1.0E-21,
    "as": 1.0E-18,
    "fs": 1.0E-15,
    "ps": 1.0E-12,
    "ns": 1.0E-09,
    "μs": 1.0E-06,
    "us": 1.0E-06,
    "ms": 1.0E-03,
    "s" : 1.0,
    "m" : 60.0,
    "h" : 60.0*60.0,
    "d" : 60.0*60.0*24.0,
    "y" : 1.0E00 * YEAR_CGS,
    "ky": 1.0E03 * YEAR_CGS,
    "My": 1.0E06 * YEAR_CGS,
    "Gy": 1.0E09 * YEAR_CGS,
    "Ty": 1.0E12 * YEAR_CGS,
    "Py": 1.0E15 * YEAR_CGS,
    "Ey": 1.0E18 * YEAR_CGS,
    "Zy": 1.0E21 * YEAR_CGS,
    "Yy": 1.0E24 * YEAR_CGS
}



def convert_half_life_readable(half_life: float) -> str:
    """
    Converts the half-life value to a human-readable format with the most suitable time unit.
    """
    if not isinstance(half_life, (int, float)) or half_life < 0:
        raise ValueError(f"Invalid input for half-life: {half_life}. It must be a positive number.")
    
    if half_life == np.inf:
        return 'Stable', np.inf, 's'

    for unit, factor in sorted(UNITS.items(), key=lambda x: x[1], reverse=True):
        if half_life >= factor:
            half_life_value = half_life / factor
            half_life_unit = unit
            return f"{half_life_value:.2f} {unit}", half_life_value, half_life_unit

    half_life_value = half_life
    half_life_unit = 's'
    return f"{half_life_value:.2f} s", half_life_value, half_life_unit




def convert_str_to_list(df: pd.DataFrame, column_name: str) -> pd.DataFrame:
    for idx, value in df[column_name].items():
        if isinstance(value, str):
            try:
                if value.startswith('[') and value.endswith(']'):
                    df.at[idx, column_name] = ast.literal_eval(value)
                else:
                    print(f"Warning: '{value}' is not a valid list string at index {idx}")
            except (ValueError, SyntaxError) as e:
                print(f"Warning: Failed to convert value at index {idx}: '{value}' ({str(e)})")
    return df





def convert_to_rational(number):
    """
    Converts a string representation of a number to a SymPy Rational object.
    """
    number = str(number).strip()

    # Handle scientific notation (e.g., 1e2, 3E-2)
    if 'e' in number or 'E' in number:
        base, exp = number.lower().split('e')
        exp = int(exp)
    else:
        base, exp = number, 0

    # Process the base part (e.g., for 123.456 or 123)
    base_parts = base.split('.')
    if len(base_parts) == 1:
        base_parts.append('')

    # Remove leading zeros in the base
    base_parts[0] = base_parts[0].lstrip('0')
    if len(base_parts[0]) == 0:
        base_parts[1] = base_parts[1].lstrip('0')

    # Combine the integer and fractional parts of the base
    numerator = base_parts[0] + base_parts[1]
    denominator = 10**len(base_parts[1])

    # Convert the number to a rational (fraction)
    rational_number = S(numerator) / S(denominator)

    # Apply the exponent (multiply or divide by 10^exp)
    if exp > 0:
        rational_number *= 10**exp
    elif exp < 0:
        rational_number /= 10**(-exp)

    return rational_number









