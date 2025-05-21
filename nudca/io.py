# -*- coding: utf-8 -*-

import re
import csv
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Any, Optional, Union, Dict, List

from .nuclide import Nuclide, Z_to_element
from .decay_database import load_decay_database


class IO:

    def __init__(self):
        pass



def add_dictionaries(dict_a: Dict[str, float], dict_b: Dict[str, float]) -> Dict[str, float]:

    new_dict = dict_a.copy()
    for nuclide, quantity in dict_b.items():
        if nuclide in new_dict:
            new_dict[nuclide] = new_dict[nuclide] + quantity
        else:
            new_dict[nuclide] = quantity

    return new_dict


def sort_dictionary_alphabetically(input_dict: Dict[str, float]) -> Dict[str, float]:

    return dict(sorted(input_dict.items(), key=lambda x: x[0]))


def _read_csv_file(filepath: Union[str, Path], delimiter: str, encoding: str) -> List[List[str]]:
    """Read CSV file. All file read side-effect is here to assist testing read_csv() with mock."""
    with open(filepath, "r", encoding=encoding) as file:
        reader_object = csv.reader(file, delimiter=delimiter)
        lines = list(reader_object)

    return lines


def _write_csv_file(filename: Union[str, Path], rows: list[list[str]], delimiter: str, encoding: str) -> None:
    """Write a CSV file from a list of rows."""
    with open(filename, "w", encoding=encoding, newline="") as file:
        writer_object = csv.writer(file, delimiter=delimiter)
        writer_object.writerows(rows)



def _parse_row(row: List[str]) -> Dict[Union[str, int], float]:
    """Parse one row that should be nuclide, quantity format."""

    if len(row) not in {2, 3}:
        raise ValueError(
            f"This row of input file ('{row}') does not satisfy required format: i.e. nuclide, "
            "quantity (, unit)."
        )

    nuc = row[0]
    nuclide: Union[str, int] = int(nuc) if nuc.isnumeric() else nuc
    quantity = float(row[1])

    return {nuclide: quantity}




def extract_mass_number(nuclide):
    match = re.search(r'\d+', nuclide)
    return int(match.group()) if match else float('inf')


def get_sorted_Y(nuclide_Y):
    sorted_Y = dict(sorted(nuclide_Y.items(), key=lambda x: extract_mass_number(x[0])))
    A, Y = [], []
    for key, value in sorted_Y.items():
        A_decay = Nuclide(key).A
        Y_decay = value
        A.append(A_decay)
        Y.append(Y_decay)
    return A, Y



def to_dict(Z, A, Y):
    nuclide_Y = {}
    for z, a, y in zip(Z, A, Y):
        element = Z_to_element(z)
        nuclide = Nuclide(f'{element}{int(a)}').nuclide_symbol
        nuclide_Y[nuclide] = y

    return nuclide_Y



def filter_Y(Z, A, Y):
    nuclide_dict = to_dict(Z, A, Y)
    decay_database = load_decay_database()
    nuclides = decay_database.nuclides
    common_nuclides = set(nuclide_dict.keys()) & set(nuclides)
    filtered_dict = {key: nuclide_dict[key] for key in common_nuclides}
    
    return filtered_dict



def group_and_sum(key_name, key_data, abundance):
    """Group data by key and sum abundances."""
    df = pd.DataFrame({key_name: key_data, 'Y': abundance})
    result = df.groupby(key_name)['Y'].sum().reset_index()
    return np.array(result[key_name]), np.array(result['Y'])


