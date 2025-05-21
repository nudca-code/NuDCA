# -*- coding: utf-8 -*-

from typing import Union


class Nuclide:
    """
    Represents a nuclear nuclide, providing properties for atomic number (Z),
    mass number (A), neutron number (N), metastable state, and canonical id.
    Supports parsing from string or integer formats and provides equality/hash operations.
    """

    def __init__(
            self,
            nuclide: Union[str, int]
            ) -> None:
        """
        Initialize a Nuclide instance from a string or integer representation.
        Args:
            nuclide (str or int): Nuclide symbol (e.g. 'U-235', '235U', 922350000) or canonical id.
        """
        self.nuclide_symbol = parse_nuclide(nuclide)


    @property
    def Z(self) -> int:
        """
        Returns the atomic number (number of protons) of the nuclide.
        Returns:
            int: Atomic number Z.
        """
        return element_to_Z(self.nuclide_symbol.split("-", maxsplit=1)[0])


    @property
    def A(self) -> int:
        """
        Returns the mass number (number of protons + neutrons) of the nuclide.
        Returns:
            int: Mass number A.
        """
        return int(self.nuclide_symbol.split("-")[1].strip('MNOPQ'))


    @property
    def N(self) -> int:
        """
        Returns the neutron number (N = A - Z) of the nuclide.
        Returns:
            int: Neutron number N.
        """
        return self.A - self.Z
        

    @property
    def state(self) -> str:
        """
        Returns the metastable state character, i.e. '' for ground state, 'M' for first metastable
        state, 'N' for second metastable state, etc.
        Returns:
            str: Metastable state symbol ('' for ground, 'M', 'N', ... for excited states).
        """
        return self.nuclide_symbol.split("-")[1].strip("0123456789")


    @property
    def id(self) -> int:
        """
        Returns the canonical nuclide id, in zzzaaassss form. Ground state is 0000, first excited
        state ("m") is 0001, second ("n") is 0002, etc.
        Returns:
            int: Canonical nuclide id.
        """
        return _build_nuclide_id(self.Z, self.A, self.state)


    def __repr__(self) -> str:
        """
        Returns a string representation of the Nuclide instance.
        Returns:
            str: Human-readable nuclide string.
        """
        rep = f"Nuclide: {self.nuclide_symbol}"
        return rep


    def __eq__(self, other: object) -> bool:
        """
        Check whether two ``Nuclide`` instances are equal with ``==`` operator.
        Args:
            other (object): Another Nuclide instance.
        Returns:
            bool: True if equal, False otherwise.
        """
        if not isinstance(other, Nuclide):
            return NotImplemented
        return self.nuclide_symbol == other.nuclide
    
    
    def __ne__(self, other: object) -> bool:
        """
        Check whether two ``Nuclide`` instances are not equal with ``!=`` operator.
        Args:
            other (object): Another Nuclide instance.
        Returns:
            bool: True if not equal, False otherwise.
        """
        if not isinstance(other, Nuclide):
            return NotImplemented
        return not self.__eq__(other)


    def __hash__(self) -> int:
        """
        Hash function for ``Nuclide`` instances, allowing use in sets and as dict keys.
        Returns:
            int: Hash value based on nuclide symbol.
        """
        return hash(self.nuclide_symbol)


class NuclideStrError(Exception):
    """
    Custom exception class for invalid nuclide strings.

    Attributes
    ----------
    nuclide_symbol : str
        The invalid nuclide string.
    error_message : str
        Contextual information about the error.
    """

    def __init__(self, nuclide: str, error_message: str, *args: object) -> None:
        """
        Initialize the exception with the invalid nuclide string and error message.
        Args:
            nuclide (str): The invalid nuclide string.
            error_message (str): Contextual information about the error.
            *args (object): Additional arguments for the base Exception class.
        """
        super().__init__(error_message, *args)
        self.nuclide_symbol = nuclide
        self.error_message = error_message


    def __str__(self) -> str:
        """
        Return a string representation of the error.
        Returns:
            str: A detailed error message.
        """
        return f"Invalid nuclide string '{self.nuclide_symbol}': {self.error_message}"


    def __repr__(self) -> str:
        """
        Developer-friendly representation of the exception.
        Returns:
            str: Debug-friendly string representation.
        """
        return (
            f"NuclideStrError(nuclide_symbol='{self.nuclide_symbol}', "
            f"error_message='{self.error_message}')"
        )
        


def _build_nuclide_string(Z: int, A: int, state: str = "") -> str:
    """
    Builds a nuclide string from given atomic mass and number.
    Args:
        Z (int): Atomic number.
        A (int): Mass number.
        state (str): Metastable state symbol (optional).
    Returns:
        str: Nuclide string in the format 'Element-A[state]'.
    """
    if Z not in ATOMIC_MASS_DICT:
        raise ValueError(f"{Z} is not a valid atomic number")
    
    return f"{ATOMIC_MASS_DICT[Z]}-{A}{state}"



def _build_nuclide_id(Z: int, A: int, state: str = "") -> int:
    """
    Builds a canonical nuclide id from atomic number, atomic mass, and energy state.
    Args:
        Z (int): Atomic number.
        A (int): Mass number.
        state (str): Metastable state symbol (optional).
    Returns:
        int: Canonical nuclide id in zzzaaassss form.
    """
    if state != "":
        if state in get_metastable_symbols():
            state_int = get_metastable_symbols().index(state) + 1
        else:
            raise ValueError(f"{state} is not a valid energy state.")
    else:
        state_int = 0
    nuclide_id = (Z * 10000000) + (A * 10000) + state_int
    return nuclide_id


def _parse_nuclide_string(nuclide_string: str) -> str:
    """
    Parses a nuclide string from e.g. '241Pu' or 'Pu241' format to 'Pu-241' format. Note this
    function works for both radioactive and stable nuclides.
    Args:
        nuclide_string (str): Input nuclide string (e.g. '241Pu', 'Pu241', 'U-235').
    Returns:
        str: Standardized nuclide string in 'Element-A[state]' format.
    Raises:
        NuclideStrError: If the string is not a valid nuclide representation.
    """
    nuclide_symbol = nuclide_string.replace(" ", "").replace("-", "", 1)  # Remove spaces and first hyphen
    # Check if the input is alphanumeric
    if not nuclide_symbol.isalnum():
        raise NuclideStrError(
            nuclide_string,
            "Nuclide strings must contain only letters, numbers, and (optionally) up to one hyphen.",
        )
    mass_number = "".join(filter(str.isdigit, nuclide_symbol))  # Extract mass number
    if len(mass_number) == 0 or int(mass_number) > 500:
        raise NuclideStrError(nuclide_string, f"Mass number ({mass_number}) is unphysical.")
    # Split the nuclide string into element symbol and metastable character
    parts = nuclide_symbol.split(mass_number)
    if len(parts) != 2:
        raise NuclideStrError(nuclide_string, "Invalid split of element and mass number.")
    element_symbol, metastable_symbol = (parts[1], parts[0]) if parts[0] == "" else (parts[0], parts[1])
    # Capitalize the element and validate it
    element_symbol = element_symbol.capitalize()
    if element_symbol not in ELEMENT_SYMBOL_DICT:
        raise NuclideStrError(nuclide_string, f"Element ({element_symbol}) is invalid.")
    # Validate metastable state character
    if len(metastable_symbol) > 1 or (metastable_symbol and metastable_symbol not in get_metastable_symbols()):
        raise NuclideStrError(nuclide_string, f"Metastable state specification ({metastable_symbol}) is invalid.")
    # Return the properly formatted nuclide string
    return f"{element_symbol}-{mass_number}{metastable_symbol}"


def _parse_nuclide_id(nuclide_id: int) -> str:
    """
    Parses a canonical nuclide id (zzzaaassss) into a nuclide string.
    Args:
        nuclide_id (int): Canonical nuclide id.
    Returns:
        str: Nuclide string in 'Element-A[state]' format.
    """
    id_zzzaaa = int(nuclide_id / 10000)
    state_digits = nuclide_id - (id_zzzaaa * 10000)
    state = get_metastable_symbols()[state_digits - 1] if state_digits > 0 else ""
    Z = int(id_zzzaaa / 1000)
    A = id_zzzaaa - (Z * 1000)
    nuclide_symbol = _build_nuclide_string(Z, A, state)
    return nuclide_symbol


def parse_nuclide(nuclide: Union[str, int]) -> str:
    """
    Parses a nuclide string or canonical id into element symbol - mass number format.
    Args:
        nuclide (str or int): Nuclide string or canonical id.
    Returns:
        str: Standardized nuclide string in 'Element-A[state]' format.
    Raises:
        TypeError: If input is not int or str.
    """
    if isinstance(nuclide, int):
        nuclide_char = _parse_nuclide_id(nuclide)
    elif isinstance(nuclide, str):
        nuclide_char = nuclide
    else:
        raise TypeError("Invalid input type, expected int or str")
    nuclide_symbol = _parse_nuclide_string(nuclide_char)
    return nuclide_symbol

         
# Atomic number to element symbol mapping
ATOMIC_MASS_DICT = {
    0  : 'Nn',
    1  : "H" ,
    2  : "He",
    3  : "Li",
    4  : "Be",
    5  : "B" ,
    6  : "C" ,
    7  : "N" ,
    8  : "O" ,
    9  : "F" ,
    10 : "Ne",
    11 : "Na",
    12 : "Mg",
    13 : "Al",
    14 : "Si",
    15 : "P" ,
    16 : "S" ,
    17 : "Cl",
    18 : "Ar",
    19 : "K" ,
    20 : "Ca",
    21 : "Sc",
    22 : "Ti",
    23 : "V" ,
    24 : "Cr",
    25 : "Mn",
    26 : "Fe",
    27 : "Co",
    28 : "Ni",
    29 : "Cu",
    30 : "Zn",
    31 : "Ga",
    32 : "Ge",
    33 : "As",
    34 : "Se",
    35 : "Br",
    36 : "Kr",
    37 : "Rb",
    38 : "Sr",
    39 : "Y" ,
    40 : "Zr",
    41 : "Nb",
    42 : "Mo",
    43 : "Tc",
    44 : "Ru",
    45 : "Rh",
    46 : "Pd",
    47 : "Ag",
    48 : "Cd",
    49 : "In",
    50 : "Sn",
    51 : "Sb",
    52 : "Te",
    53 : "I" ,
    54 : "Xe",
    55 : "Cs",
    56 : "Ba",
    57 : "La",
    58 : "Ce",
    59 : "Pr",
    60 : "Nd",
    61 : "Pm",
    62 : "Sm",
    63 : "Eu",
    64 : "Gd",
    65 : "Tb",
    66 : "Dy",
    67 : "Ho",
    68 : "Er",
    69 : "Tm",
    70 : "Yb",
    71 : "Lu",
    72 : "Hf",
    73 : "Ta",
    74 : "W" ,
    75 : "Re",
    76 : "Os",
    77 : "Ir",
    78 : "Pt",
    79 : "Au",
    80 : "Hg",
    81 : "Tl",
    82 : "Pb",
    83 : "Bi",
    84 : "Po",
    85 : "At",
    86 : "Rn",
    87 : "Fr",
    88 : "Ra",
    89 : "Ac",
    90 : "Th",
    91 : "Pa",
    92 : "U" ,
    93 : "Np",
    94 : "Pu",
    95 : "Am",
    96 : "Cm",
    97 : "Bk",
    98 : "Cf",
    99 : "Es",
    100: "Fm",
    101: "Md",
    102: "No",
    103: "Lr",
    104: "Rf",
    105: "Db",
    106: "Sg",
    107: "Bh",
    108: "Hs",
    109: "Mt",
    110: "Ds",
    111: "Rg",
    112: "Cn",
    113: "Nh",
    114: "Fl",
    115: "Mc",
    116: "Lv",
    117: "Ts",
    118: "Og",
}
# Reverse mapping: element symbol to atomic number
ELEMENT_SYMBOL_DICT = dict((value, key) for key, value in ATOMIC_MASS_DICT.items())
# List of allowed metastable state symbols
METASTABLE_CHARS = ['M', 'N', 'O', 'P', 'Q', 'R']


def get_metastable_symbols() -> list[str]:
    """
    Returns the list of allowed metastable state symbols.
    Returns:
        list[str]: List of metastable state characters.
    """
    return METASTABLE_CHARS


def Z_to_element(Z: int) -> str:
    """
    Convert atomic number to element symbol.
    Args:
        Z (int): Atomic number.
    Returns:
        str: Element symbol.
    Raises:
        ValueError: If atomic number is invalid.
    """
    try:
        return ATOMIC_MASS_DICT[Z]
    except KeyError:
        raise ValueError(f"Invalid atomic number: {Z}")


def element_to_Z(element_symbol: str) -> int:
    """
    Convert element symbol to atomic number.
    Args:
        element_symbol (str): Element symbol (e.g. 'U', 'Pu').
    Returns:
        int: Atomic number.
    Raises:
        ValueError: If element symbol is invalid.
    """
    try:
        return ELEMENT_SYMBOL_DICT[element_symbol]
    except KeyError:
        raise ValueError(f"Invalid element symbol: {element_symbol}")
