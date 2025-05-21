
import os
import sys
from pathlib import Path, PurePath
ENDFtk_path = Path(__file__).resolve().parents[1].joinpath('ENDFtk', 'build')
sys.path.insert(0, str(ENDFtk_path))
sys.path.insert(1, str(Path(__file__).resolve().parents[1]))

import re
import ast
import logging
import numpy as np
import pandas as pd
from sympy import Rational

import zipfile
import requests
from bs4 import BeautifulSoup
from concurrent.futures import ThreadPoolExecutor


import ENDFtk
from nudca.nuclide import Nuclide, Z_to_element
from nudca.constants import N_MASS, ATOMIC_MASS, YEAR_CGS
# from src.decay.utils.converters import convert_half_life_readable


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



def get_nuclide_symbol(filename):
    """Extract nuclide symbol from filename."""
    match = re.search(r'([A-Za-z]+)-(\d+)([A-Z]?)', filename)
    if match:
        element_symbol, mass_number, metastable = match.groups()
        return Nuclide(f'{element_symbol}-{mass_number}{metastable}').nuclide_symbol
    else:
        raise ValueError("Invalid filename format. Expected 'Element-Mass[Metastable]'.")


#-------------------------------------------------------------------#
#                           ENDFDownloader
#-------------------------------------------------------------------#

class ENDFDownloader:
    def __init__(self, base_url, download_folder="downloads"):
        self.base_url = base_url
        self.download_folder = download_folder

        os.makedirs(self.download_folder, exist_ok=True)

        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s - %(levelname)s - %(message)s"
        )

    def fetch_webpage(self, url):
        try:
            response = requests.get(url, timeout=10)
            response.raise_for_status()
            return response.text
        except requests.RequestException as e:
            logging.error(f"Error fetching the URL {url}: {e}")
            return None


    def download_task(self, file_url, save_path):
        try:
            response = requests.get(file_url, stream=True, timeout=20)
            response.raise_for_status()
            with open(save_path, "wb") as file:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        file.write(chunk)
            logging.info(f"Downloaded: {file_url}")
        except requests.RequestException as e:
            logging.error(f"Error downloading the file {file_url}: {e}")


    def unzip_file(self, zip_path):
        try:
            with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                for file in zip_ref.namelist():
                    base_name, ext = os.path.splitext(file)
                    extracted_path = os.path.join(self.download_folder, base_name + ".endf")
                    counter = 1
                    while os.path.exists(extracted_path):
                        extracted_path = os.path.join(
                            self.download_folder, f"{base_name}_{counter}.endf"
                        )
                        counter += 1
                    with open(extracted_path, "wb") as f:
                        f.write(zip_ref.read(file))
                    logging.info(f"Extracted and renamed: {extracted_path}")
        except zipfile.BadZipFile as e:
            logging.error(f"Error unzipping the file {zip_path}: {e}")


    def process_file(self, href):
        file_url = self.base_url + href
        zip_save_path = os.path.join(self.download_folder, href)
        self.download_task(file_url, zip_save_path)
        self.unzip_file(zip_save_path)

        try:
            os.remove(zip_save_path)
            logging.info(f"Deleted ZIP file: {zip_save_path}")
        except OSError as e:
            logging.error(f"Error deleting ZIP file {zip_save_path}: {e}")


    def download_files(self):
        html_content = self.fetch_webpage(self.base_url)
        if not html_content:
            return

        soup = BeautifulSoup(html_content, "html.parser")
        tasks = [link.get("href") for link in soup.find_all("a", href=True) if link.get("href").endswith(".zip")]

        with ThreadPoolExecutor(max_workers=int(3/4 * os.cpu_count())) as executor:
            futures = [executor.submit(self.process_file, task) for task in tasks]
            for future in futures:
                try:
                    future.result()
                except Exception as e:
                    logging.error(f"Error processing file: {e}")




#-------------------------------------------------------------------#
#                           ENDFFileHandle
#-------------------------------------------------------------------#

class ENDFFileHandle:
    def __init__(self, file_path):
        if not Path(file_path).is_file():
            raise FileNotFoundError(f"File not found: {file_path}")
        self._file_path = file_path


    def _load_tape(self):
        """Load the ENDF file into a Tape object."""
        return ENDFtk.tree.Tape.from_file(self._file_path)


    def _get_material_section(self, file_number, section_number):
        """Retrieve a specific section from the material."""
        tape = self._load_tape()
        material_number = self.get_material_number()
        return tape.material(material_number).file(file_number).section(section_number).parse()


    def print_MAT_MF(self):
        """Print the MAT and MF structure."""
        tape = self._load_tape()
        for material in tape.materials:
            print(f"Material {material.material_number:4}")
            for file in material.files:
                print(f"  File {file.file_number:2}")
                for section in file.sections:
                    print(f"    Section {section.section_number:3}, {section.NC:5} lines")


    def get_material_number(self):
        """Get the material number of the file."""
        tape = self._load_tape()
        for material in tape.materials:
            material_number = material.material_number
        return material_number


    def get_file_numbers(self):
        """Get all file numbers in the material."""
        tape = self._load_tape()
        return [file.file_number for material in tape.materials for file in material.files]


    def print_file1_description(self):
        """Print the description from file 1, section 451."""
        section = self._get_material_section(1, 451)
        print(section.description)


    def get_Z_and_A(self):
        """Get Z and A (atomic number and mass number)."""
        section = self._get_material_section(1, 451)
        ZA = section.ZA
        Z, A = divmod(ZA, 1000)
        return Z, A

    def get_half_life(self):
        """Get the half-life of the nuclide."""
        section_451 = self._get_material_section(1, 451)
        if section_451.NSUB == 4:
            section_457 = self._get_material_section(8, 457)
            half_life = section_457.half_life
            # return getattr(section_457, 'half_life', None)
        return np.float64(half_life)
    

    def get_atomic_weight_ratio(self):
        section_451 = self._get_material_section(1, 451)
        if section_451.NSUB == 4:
            section_457 = self._get_material_section(8, 457)
            return getattr(section_457, 'atomic_weight_ratio', None)  # AWR
        return None



    def is_stable(self):
        """Check if the nuclide is stable."""
        section = self._get_material_section(8, 457)
        return getattr(section, 'is_stable', False)


    def get_num_decay_modes(self):
        """Get the number of decay modes."""
        section = self._get_material_section(8, 457)
        return section.decay_modes.number_decay_modes


    @staticmethod
    def decay_mode_type(RTYP):
        flags = {
            0:      'γ',
            1:      'β-',
            2:      'β+&EC',
            3:      'IT',
            4:      'α',
            5:      'n',
            6:      'SF',
            7:      'p',
            10:     'Unknown',
            1.1:    '2β-',
            1.4:    'β-,α',
            1.5:    'β-,n',
            1.55:   'β-,2n',
            1.555:  'β-,3n',
            1.5555: 'β-,4n',
            2.4:    'β+,α',
            2.5:    'β+,n',
            2.6:    'β+,SF',
            2.7:    'β+,p',
            2.77:   'β+,2p',
            5.5:    '2n',
            7.7:    '2p'
        }
        return flags.get(RTYP, 'Other')


    def decay_radiation_type(STYP):
        flags = {
            0:  'γ',
            1:  'β-',
            2:  'β+&EC',
            4:  'α',
            5:  'n',
            6:  'SF',
            7:  'p',
            8:  'e-',
            9:  'x',
            10: 'anue',
            11: 'nue',
        }
        return flags.get(STYP, 'Other')



    def get_decay_modes(self):
        """Get all decay modes."""
        section = self._get_material_section(8, 457)

        # return [mode.RTYP for mode in section.decay_modes.decay_modes]
        return [self.decay_mode_type(mode.RTYP) for mode in section.decay_modes.decay_modes]


    def get_decay_progeny(self):
        Z, A = self.get_Z_and_A()
        section = self._get_material_section(8, 457)

        decay_type_mapping = {
            0:      (Z, A),             # γ
            1:      (Z + 1, A),         # β-
            2:      (Z - 1, A),         # β+&EC
            3:      (Z, A),             # IT
            4:      (Z - 2, A - 4),     # α
            5:      (Z, A - 1),         # n
            7:      (Z - 1, A - 1),     # p
            1.1:    (Z + 2, A),         # 2β-
            1.4:    (Z - 1, A - 4),     # β-α
            1.5:    (Z + 1, A - 1),     # β-n
            1.55:   (Z + 1, A - 2),     # β-,2n
            1.555:  (Z + 1, A - 3),     # β-,3n
            1.5555: (Z + 1, A - 4),     # β-,4n
            2.4:    (Z - 3, A - 4),     # β+,α
            2.5:    (Z - 1, A - 1),     # β+,n
            2.7:    (Z - 2, A - 1),     # β+,p
            2.77:   (Z - 3, A - 2),     # β+,2p
            5.5:    (Z, A - 2),         # 2n
            7.7:    (Z - 2, A - 2),     # 2p
        }

        nuclides_progeny = []
        for mode in section.decay_modes.decay_modes:
            decay_type = mode.RTYP
            state_daughter = mode.RFS
            
            if decay_type in decay_type_mapping:
                Z_progeny, A_progeny = decay_type_mapping[decay_type]
                element_symbol = Z_to_element(Z_progeny)
                metastable = self.metastable_char(state_daughter)
                nuclide_progeny = f"{element_symbol}-{A_progeny}{metastable}"
            else:
                nuclide_progeny = 'X'

            nuclides_progeny.append(nuclide_progeny)

        return nuclides_progeny
        

    def get_branching_ratios(self):
        """Get branching ratios for decay modes."""
        section = self._get_material_section(8, 457)
        return [np.float64(mode.BR[0]) for mode in section.decay_modes.decay_modes]


    def get_decay_chains(self):
        """Get decay chains."""
        section = self._get_material_section(8, 457)
        return [mode.decay_chain for mode in section.decay_modes.decay_modes]


    def get_parent_isomeric_states(self):
        """Get parent isomeric states."""
        section = self._get_material_section(8, 457)
        return section.LISO
    

    @staticmethod
    def metastable_char(RFS):
        state = RFS
        if state == 0.0:
            metastable_char = ''
        elif state == 1.0:
            metastable_char = 'M'
        elif state == 2.0:
            metastable_char = 'N'
        elif state == 3.0:
            metastable_char = 'O'
        elif state == 4.0:
            metastable_char = 'P'
        elif state == 4.0:
            metastable_char = 'Q'
        else:
            raise ValueError('Be Careful!')
        return metastable_char
    

    def get_parent_metastable_state(self):
        section = self._get_material_section(8, 457)
        return self.metastable_char(section.LISO)


    def get_daughter_isomeric_states(self):
        """Get daughter isomeric states."""
        section = self._get_material_section(8, 457)
        return [mode.RFS for mode in section.decay_modes.decay_modes]
    

    def get_decay_q_values(self):
        """Get Q values for decay modes."""
        section = self._get_material_section(8, 457)
        return [mode.q_value[0] for mode in section.decay_modes.decay_modes]


    def get_effective_Q(self):
        """Get effective Q value for decay modes."""
        section = self._get_material_section(1, 451)

        pattern = r"Q effective:\s+([\d.E+-]+)\s+keV"
        match = re.search(pattern, section.description)
        
        if match:
            Q_effective = float(match.group(1)) * 1.e3  # Convert to eV
            return Q_effective
        else:
            return 0



    def get_mean_energy(self, energy_marker):
        section = self._get_material_section(1, 451)

        pattern = f'{energy_marker}' + r":\s+([\d.E+-]+)\s+\+-\s+([\d.E+-]+)\s+keV"
        match = re.search(pattern, section.description)
        
        if match:
            mean_energy = float(match.group(1)) * 1.e3  # Convert to eV
            uncertainty = float(match.group(2)) * 1.e3 
            return mean_energy
        else:
            return 0.0
        

    def get_mean_neutrino_energy(self):
        """Get neutrino energy for decay modes."""
        return self.get_mean_energy('Mean Neutrino Energy')
        

    def get_mean_gamma_energy(self):
        return self.get_mean_energy('Mean Gamma Energy')
    

    def get_mean_xray_energy(self):
        return self.get_mean_energy('Mean X-Ray+511 Energy')
    

    def get_mean_Bminus_energy(self):
        return self.get_mean_energy('Mean B- Energy')
    

    def get_mean_Bplus_energy(self):
        return self.get_mean_energy('Mean B+ Energy')
    

    def get_mean_alpha_energy(self):
        return self.get_mean_energy('Mean Alpha Energy')
    

    def get_mean_neutron_energy(self):
        return self.get_mean_energy('Mean Neutron Energy')
    

    def get_mean_proton_energy(self):
        return self.get_mean_energy('Mean Proton Energy')


    
    def get_number_decay_energies(self):
        """Get the number of decay modes."""
        section = self._get_material_section(8, 457)
        return section.average_decay_energies.number_decay_energies
    

    def get_electromagnetic_decay_energy(self):
        section = self._get_material_section(8, 457)
        num_E = section.average_decay_energies.number_decay_energies
        if num_E == 3:
            decay_energy = section.average_decay_energies.electromagnetic_decay_energy[0]

        return np.float64(decay_energy)

    
    def get_heavy_particle_decay_energy(self):
        section = self._get_material_section(8, 457)
        num_E = section.average_decay_energies.number_decay_energies
        if num_E == 3:
            decay_energy = section.average_decay_energies.heavy_particle_decay_energy[0]

        return np.float64(decay_energy)
    

    def get_light_particle_decay_energy(self):
        section = self._get_material_section(8, 457)
        num_E = section.average_decay_energies.number_decay_energies
        if num_E == 3:
            decay_energy = section.average_decay_energies.light_particle_decay_energy[0]

        return np.float64(decay_energy)



    def test(self):
        tape = ENDFtk.tree.Tape.from_file(self._file_path)
        material_number = self.get_material_number()
        xs = tape.material(material_number).file(1).section(451).parse()
        # xs = tape.material(material_number).file(8).section(457).parse()

        # print(dir(xs))

        # print(xs.description)

        print(self.get_effective_Q())

        print(self.get_parent_isomeric_states())


        # print(dir(xs.average_decay_energies.E))

        # print(xs.average_decay_energies.number_decay_energies)
            


        # print(list(xs.average_decay_energies.E)
              
        # for energy in xs.average_decay_energies.E:
        #     print(list(energy))
        

        # description = self.print_file1_description()
        # print(dir(xs.average_decay_energies))
        # print(xs.average_decay_energies.number_decay_energies)

        # print(xs.average_decay_energies.heavy_particle_decay_energy[0])

        # for energy in xs.average_decay_energies.decay_energies:
        #     print(list(energy))

        # print(xs.decay_spectra)

        # for spectra in xs.decay_spectra:
        #     print(dir(spectra))

        #     print(spectra.NC, spectra.STYP, spectra.radiation_type)

            # print(spectra.average_decay_energy)



#-----------------------------------------------------------------------#
#
#-----------------------------------------------------------------------#

def process_endf_file(filepath):
    """
    Process a single ENDF file and return the data to be added to the DataFrame.
    """
    endftk = ENDFFileHandle(filepath)
    Z, A = endftk.get_Z_and_A()

    atomic_mass_ratio = endftk.get_atomic_weight_ratio()
    atomic_mass = atomic_mass_ratio * N_MASS / ATOMIC_MASS
    half_life_s, half_life_err = endftk.get_half_life()
    half_life_readable, half_life_value, half_life_unit = convert_half_life_readable(half_life_s)
    
    decay_chains = endftk.get_decay_chains()


    return {
        'Radionuclide': get_nuclide_symbol(filepath),
        'Element': Z_to_element(Z),
        'Z': Z,
        'A': A,
        'Metastable_State': endftk.get_parent_metastable_state(),
        'Atomic_Mass': atomic_mass,
        'Is_Stable': endftk.is_stable(),
        'Half_Life_Second': half_life_s,
        'Half_Life_Readable': half_life_readable,
        'Half_Life_Value': half_life_value,
        'Half_Life_Unit': half_life_unit,
        # 'Half_Life_Err': half_life_err,
        'Num_Decay_Modes': endftk.get_num_decay_modes(),
        'Decay_Modes': endftk.get_decay_modes(),
        'Branching_Ratios': endftk.get_branching_ratios(),
        'Progeny': endftk.get_decay_progeny(),
        'Decay_Energy_EM': endftk.get_electromagnetic_decay_energy(),
        'Decay_Energy_LP': endftk.get_light_particle_decay_energy(),
        'Decay_Energy_HP': endftk.get_heavy_particle_decay_energy(),
        'Decay_Energy_Neutrino': endftk.get_mean_neutrino_energy(),
        'Decay_Energy_Gamma': endftk.get_mean_gamma_energy(),
        'Decay_Energy_Xray': endftk.get_mean_xray_energy(),
        'Decay_Energy_Beta_Minus': endftk.get_mean_Bminus_energy(),
        'Decay_Energy_Beta_Plus': endftk.get_mean_Bplus_energy(),
        'Decay_Energy_Alpha': endftk.get_mean_alpha_energy(),
        'Decay_Energy_Neutron': endftk.get_mean_neutron_energy(),
        'Decay_Energy_Proton': endftk.get_mean_proton_energy(),
        'Decay_Effective_Q': endftk.get_effective_Q(),
        # 'Decay_Q_Values': endftk.get_decay_q_values(),
        # 'Isomeric_States': endftk.get_parent_isomeric_states(),
        # 'Daughter_Isomeric_States': endftk.get_daughter_isomeric_states(),
    }



def gen_endf_database(data_dirpath, database_name = 'endf_test'):
    endf_files = [f for f in os.listdir(data_dirpath) if os.path.isfile(os.path.join(data_dirpath, f))]
    
    # Process all files and accumulate the data in a list
    data = []
    for filename in endf_files:
        filepath = os.path.join(data_dirpath, filename)
        data.append(process_endf_file(filepath))
        
    df = pd.DataFrame(data)

    # Normalize Atomic_Mass using C-12 atomic mass
    mass_c12 = df[(df['Z'] == 6) & (df['A'] == 12)]['Atomic_Mass'].values
    if mass_c12.size == 1:
        ratio = 12.0 / mass_c12[0]
        df['Atomic_Mass'] *= ratio
    else:
        raise ValueError('Invalid data for C-12 atomic mass')
    
    # df['Atomic_Mass_sympy'] = df['Atomic_Mass'].apply(lambda x: Rational(x).limit_denominator())

    df.loc[df['Half_Life_Second'] == 0, 'Is_Stable'] = True
    df.loc[df['Is_Stable'] == True, 'Half_Life_Second'] = np.inf
    df.loc[df['Is_Stable'] == True, 'Half_Life_Value'] = np.inf
    df.loc[df['Is_Stable'] == True, 'Half_Life_Readable'] = 'Stable'
    df['Progeny'] = df.apply(lambda row: [None] if row['Is_Stable'] else row['Progeny'], axis=1)
    df['Decay_Modes'] = df.apply(lambda row: [None] if row['Is_Stable'] else row['Decay_Modes'], axis=1)
    

    # df = df.sort_values(by=['A', 'Z'])
    # df.sort_values(by=['A', 'Z', 'Metastable_State'], inplace=True, ascending=[False, False, False])

    df.to_csv(f'data/{database_name}.csv', index=False)
    df.to_json(f'data/{database_name}.json', orient="records", indent=4)




if __name__ == '__main__':
    # Download ENDF file
    # base_url = "https://www-nds.iaea.org/public/download-endf/ENDF-B-VIII.1/decay/"
    # ENDFDownloader(base_url).download_files()


    # ENDF FileHandle
    endf_dirpath = Path(__file__).resolve().parents[1].joinpath('data', 'ENDF', 'ENDF-B-VIII.1')
    gen_endf_database(endf_dirpath, 'ENDF-B-VIII.1_decay')

    # endf_dirpath = Path(__file__).resolve().parents[1].joinpath('data', 'ENDF', 'ENDF-B-VII.1')
    # gen_endf_database(endf_dirpath, 'ENDF-B-VII.1_decay')


    # ENDF test
    # endf_dirpath = Path(__file__).resolve().parents[1].joinpath('data', 'ENDF', 'ENDF-B-VIII.1')
    # endf_files = [f for f in os.listdir(endf_dirpath) if os.path.isfile(os.path.join(endf_dirpath, f))]
    # for endf_file in endf_files:
    #     filepath = os.path.join(endf_dirpath, endf_file)
    #     endf = ENDFFileHandle(str(filepath))
    #     endf.test()


    # filepath = Path(__file__).resolve().parents[1].joinpath('data', 'ENDF', 'ENDF-B-VIII.1', 'decay_047-Ag-95O_1335.endf')
    # endf = ENDFFileHandle(str(filepath))
    # endf.test()
    
    # print(endf.is_stable())

    # endf.print_MAT_MF()

    # endf.print_file1_description()

    
    # print(endf.get_parent_metastable_state())
    # endf.get_decay_q_values()

    # print(endf.get_decay_modes())

    # print(endf.get_half_life())


    



