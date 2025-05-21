
import numpy as np
import pandas as pd
import fortranformat as ff


NUBASE_FORMAT = ff.FortranRecordReader('(a3,a1,a4,a3,a5,a1,a1,a13,a11,a12,a11,a2,a1,a1,a9,a2)')

class Nubase:
    
    def __init__(self, file_nubase):
        self.file = file_nubase


    def read_nubase(self):
        with open(self.file, 'r') as nubase:
            nubase_lines = nubase.readlines()

        for i in range(len(nubase_lines) -26):
            # starting from the first isotope on line 38, read each line seperately
            isotope = NUBASE_FORMAT.read(nubase_lines[i+26])
            # location in list  quantity        description
            # isotope[0]        AAA             Mass number
            #        [2]        ZZZi            Atomic number, i indicated isomer number
            #        [4]        ZZZAA           Nuclide name, mass number followed by element
            #        [5]        A               m,n (isomers); p,q (levels); r (resonance); i,j (IAS)
            #        [7]        ZZZZZZ.ZZZZZZZ  Mass excess in keV
            #        [8]        ZZZZ.ZZZZZZ     Mass excess uncertainty
            #        [9]        ZZZZZ.ZZZZZZ    Isomer excitation energy in keV
            #        [10]       ZZZZ.ZZZZZZ     Isomer excitation energy uncertainty
            #        [14]       ZZZZ.ZZZZ       Half life
            #        [15]       AA              Half life units
            #        [16]       ZZZ.ZZZ         Half life uncertainty

            meta_state = isotope[5].strip()
            # name is combination of element and number of nucleons
            element = isotope[4].strip(" 1234567890")
            A = isotope[0].lstrip("0")
            name = (element + "-" + A + meta_state)
            # split by '#' to remove characters prior to float conversion
            # convert isomer excitation energy to amu if nonzero
            isomer_excitation_energy = (
                float((isotope[9].strip()).split("#")[0]) / amu
                if isotope[9].strip(" non-exist") != ""
                else 0
                )
            

            nubase_df = pd.DataFrame()
            isotope_dictionary_list = {}
            stable_dictionary_list = {}

            # convert half-life to seconds
            half_life_raw = isotope[14].strip(" ~<>#").strip()
            # only include nuclides with half-life information
            if half_life_raw not in ["", "stbl", "p-unst"]:
                half_life = float(half_life_raw)
                half_life_units = isotope[15].strip()
                half_life_s = convert_half_life(half_life, half_life_units)
                # build dictionary for specific isotope
                isotope_dictionary_list[name] = [element, int(A), isomer_excitation_energy,
                                                half_life, half_life_units, half_life_s]
            elif name == "Ta-180m":
                # only stable isomer (no half life data, leave blank)
                isotope_dictionary_list["Ta-180m"] = ["Ta", 180, isomer_excitation_energy,
                                                    np.nan, np.nan, np.nan]
            elif half_life_raw == "stbl":
                stable_dictionary_list[name] = [element, int(A), isomer_excitation_energy,
                                                    np.nan, np.nan, np.nan]


        # convert from dictionaries to DataFrame for max speed
        nubase_df = pd.DataFrame.from_dict(isotope_dictionary_list, orient='index')
        nubase_stable_df = pd.DataFrame.from_dict(stable_dictionary_list, orient='index')
        nubase_df.columns = ["Element", "A", "Isomer_excitation_energy", "Half_life",
                            "Half_life_units", "Half_life_s"]
        nubase_stable_df.columns = ["Element", "A", "Isomer_excitation_energy", "Half_life",
                            "Half_life_units", "Half_life_s"]