
import pandas as pd
import fortranformat as ff


AME2020_FORMAT = ff.FortranRecordReader('(a1,a3,a5,a5,a5,1x,a3,a4,1x,a14,a12,a13,1x,a10,1x,a2,a13,a11,1x,a3,1x,a13,a12)')

class AME:

    def __init__(self, file_ame):
        self.file = file_ame

    
    def read_AME(self):
        with open(self.file, 'r') as ame2020:
            ame2020_lines = ame2020.readlines()

        ame2020_df = pd.DataFrame()
        isotope_dictionary_list = {}

        for i in range(len(ame2020_lines)-37):
            # starting from the first isotope on line 38, read each line seperately
            isotope = AME2020_FORMAT.read(ame2020_lines[i+37])
            # name is combination of element and number of nucleons
            name = isotope[5].split()[0] + "-" + str(isotope[4]).split()[0]
            Z = int(isotope[3].split()[0])
            A = int(isotope[4].split()[0])
            # concatenate atomic number with decimal places, in μ amu
            # split by '#' to remove characters prior to float conversion
            # divide by 10e6, as AME2020 is in micro amu
            mass = (float(isotope[14].strip("# ")
                + isotope[15].strip("# ")))/10**6
            mass_str = isotope[14].strip("# ") + "." + (isotope[15].strip("# ")).replace(".", "")
            # build dictionary for specific isotope
            isotope_dictionary_list[name] = [Z, A, mass, mass_str]

        # convert from dictionaries to DataFrame for max speed
        ame2020_df = pd.DataFrame.from_dict(isotope_dictionary_list,
                                                orient='index')
        ame2020_df.columns = ["Z", "A", "Mass", "Mass_str"]