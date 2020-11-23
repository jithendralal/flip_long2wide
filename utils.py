import json
import os
from pathlib import Path
from tkinter import filedialog


CONFIG_FILENAME = "config.json"
CONFIG_PATH = os.getcwd()
CONFIG_FILE = os.path.join(CONFIG_PATH, CONFIG_FILENAME)


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


def open_dir(current_dir):
    return filedialog.askdirectory(initialdir=current_dir, title="Please select a directory")


def get_home():
    return str(Path.home())


def get_config():
    config_json = None
    if os.path.isfile(CONFIG_FILE):
        with open(CONFIG_FILE) as f:
            config_json = json.load(f)
    else:
        with open(CONFIG_FILE, "w+") as f:
            config_json = {"cwd": get_home()}
            json.dump(config_json, f, indent=4)
    return config_json


def set_config(key, value):
    with open(CONFIG_FILE) as f:
        config = json.load(f)
    with open(CONFIG_FILE, "w") as f:
        config[key] = value
        json.dump(config, f, indent=4)
