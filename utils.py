import json
import os
from pathlib import Path
from tkinter import filedialog

BRUKER_VARIABLES = [
    'quantity_units', 
    'analyte_name', 
    'data_set', 
    'sample_type',
    'rt_min',
    'm_z_expected',
    'area_of_pi'
]

WATERS_VARIABLES = [
    'conc', 
    'analyte_name', 
    'sample_text', 
    'type',
    'rt',
    'area'
]

WATERS_HELP_VARIABLES = [
    'conc', 
    'sample_text', 
    'type',
    'rt',
    'area'
]


def get_files(dir_name, file_type):
    list_of_files = sorted(os.listdir(dir_name))
    all_files = list()
    for file in list_of_files:
        if not file.startswith('Flipped'):
            full_path = os.path.join(dir_name, file)
            if os.path.isdir(full_path):
                all_files = all_files + get_files(full_path, file_type)
            else:
                if file.endswith(file_type):
                    all_files.append(full_path)
    return all_files


def get_home():
    return str(Path.home())


def double_it(x):
    return 2 * x


def unit_conc(x, mol_weight):
    # divide concentration by the molecular weight for each analyte
    return x / mol_weight


def unit_conc_micro(x, mol_weight):
    # divide concentration by the molecular weight for each analyte
    return x / mol_weight / 1000.0
