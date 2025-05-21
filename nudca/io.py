# -*- coding: utf-8 -*-

import re
import csv
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Any, Optional, Union, Dict, List

from .nuclide import Nuclide, Z_to_element
from .decay_database import load_decay_database


class Outputer:
    """
    Output utility for saving calculation results (abundances, heating rates) to file.
    Supports CSV, Excel, and JSON formats.
    """
    def __init__(self):
        pass


    @staticmethod
    def save_abundance(
        abundance: Union[Dict[str, float], pd.DataFrame],
        filename: Union[str, Path],
        filetype: Optional[str] = None,
        index: bool = True
    ):
        """
        Save nuclide abundance (dict or DataFrame) to file.
        Args:
            abundance: Dict or DataFrame of nuclide abundances (can be time series DataFrame).
            filename: Output file path.
            filetype: 'csv', 'xlsx', or 'json'. If None, inferred from filename.
            index: Whether to write DataFrame index (default True).
        """
        if filetype is None:
            filetype = str(filename).split('.')[-1].lower()
        if isinstance(abundance, dict):
            df = pd.DataFrame(list(abundance.items()), columns=["nuclide", "abundance"])
        else:
            df = abundance
        if filetype == 'csv':
            df.to_csv(filename, index=index)
        elif filetype in ('xls', 'xlsx'):
            df.to_excel(filename, index=index)
        elif filetype == 'json':
            df.to_json(filename, orient='records', indent=2)
        else:
            raise ValueError(f"Unsupported filetype: {filetype}")

    @staticmethod
    def save_heating_rate(
        heating_rate: Union[List[float], np.ndarray, pd.Series, pd.DataFrame],
        times: Optional[Union[List[float], np.ndarray]] = None,
        filename: Union[str, Path] = "heating_rate.csv",
        filetype: Optional[str] = None
    ):
        """
        Save heating rate (with optional time axis) to file.
        Args:
            heating_rate: List/array/Series/DataFrame of heating rates.
            times: Optional time points (same length as heating_rate).
            filename: Output file path.
            filetype: 'csv', 'xlsx', or 'json'. If None, inferred from filename.
        """
        if filetype is None:
            filetype = str(filename).split('.')[-1].lower()
        if isinstance(heating_rate, (list, np.ndarray, pd.Series)):
            if times is not None:
                df = pd.DataFrame({"time": times, "heating_rate": heating_rate})
            else:
                df = pd.DataFrame({"heating_rate": heating_rate})
        else:
            df = heating_rate
        if filetype == 'csv':
            df.to_csv(filename, index=False)
        elif filetype in ('xls', 'xlsx'):
            df.to_excel(filename, index=False)
        elif filetype == 'json':
            df.to_json(filename, orient='records', indent=2)
        else:
            raise ValueError(f"Unsupported filetype: {filetype}")


class Inputer:
    """
    Input utility for reading initial nuclide abundances from file.
    Supports CSV, TSV, Excel, and JSON formats.
    """
    def __init__(self):
        pass

    @staticmethod
    def read_abundance(
        filename: Union[str, Path],
        filetype: Optional[str] = None,
        nuclide_col: str = "nuclide",
        abundance_col: str = "abundance"
    ) -> Dict[str, float]:
        """
        Read initial nuclide abundances from file.
        Args:
            filename: Input file path.
            filetype: 'csv', 'tsv', 'xlsx', or 'json'. If None, inferred from filename.
            nuclide_col: Column name for nuclide symbol (default 'nuclide').
            abundance_col: Column name for abundance (default 'abundance').
        Returns:
            Dict[str, float]: Mapping from nuclide symbol to abundance.
        """
        if filetype is None:
            filetype = str(filename).split('.')[-1].lower()
        if filetype == 'csv':
            df = pd.read_csv(filename)
        elif filetype == 'tsv':
            df = pd.read_csv(filename, sep='\t')
        elif filetype in ('xls', 'xlsx'):
            df = pd.read_excel(filename)
        elif filetype == 'json':
            df = pd.read_json(filename)
        else:
            raise ValueError(f"Unsupported filetype: {filetype}")
        # Try to find the correct columns
        cols = [col.lower() for col in df.columns]
        ncol = nuclide_col if nuclide_col in df.columns else None
        acol = abundance_col if abundance_col in df.columns else None
        if ncol is None:
            for c in df.columns:
                if c.lower() in ["nuclide", "isotope", "nuc", "species"]:
                    ncol = c
                    break
        if acol is None:
            for c in df.columns:
                if c.lower() in ["abundance", "quantity", "amount", "y"]:
                    acol = c
                    break
        if ncol is None or acol is None:
            raise ValueError(f"Could not find nuclide or abundance columns in {filename}")
        # Standardize nuclide symbols
        nuclide_abundance = {}
        for _, row in df.iterrows():
            symbol = Nuclide(row[ncol]).nuclide_symbol
            nuclide_abundance[symbol] = float(row[acol])
        return nuclide_abundance


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


