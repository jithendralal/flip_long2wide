# coding: utf-8

import json
import logging
import numpy as np
import os
import pandas as pd
import re
import sys
import janitor

import tkinter as tk
from datetime import datetime
from IPython.display import clear_output
from pathlib import Path
from pandas import ExcelWriter
from tkinter import filedialog
from tkinter import scrolledtext


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


class DataFile:
    @classmethod
    def read(cls, file, file_type):
        f = "read_" + file_type
        # for invalid file type, fn is the lambda function that returns error message
        fn = getattr(cls, f, lambda: "Invalid file type")
        return fn(file)

    def read_csv(file):
        return pd.read_csv(file)

    def read_xlsx(file):
        return pd.read_excel(file)

    def read_txt(file):
        # return pd.read_fwf(file, delimiter="\t")
        return pd.read_csv(file, sep='\t', lineterminator='\r')


class Application(tk.Frame):
    def __init__(self, master=None):
        tk.Frame.__init__(self, master)
        self.selected_cwd = False
        self.pack(fill='both', padx=2, pady=2)
        self.master.title('ANPC - Flippy')
        self.file_type = tk.StringVar(value=".xlsx")
        self.machine_type = tk.StringVar(value="waters")
        self.load_config()
        self.cwd_label_text = ""
        self.create_widgets()
        self.df = pd.DataFrame()

    def load_config(self):
        self.config = get_config()

    def get_current_dir(self):
        return self.config["cwd"]

    def fill_analyte_name(self, df):
        analytes = df[df['analyte_name'].apply(lambda x: str(x).startswith('Compound'))].index.tolist()

        print(analytes)

        analytes_index = 0
        fill_value = df['analyte_name'][analytes[analytes_index]].split(':')[1].strip()
        for i in range(2, len(df.index)):
            if not str(df['analyte_name'][i]).startswith('Compound'):
                df.at[i, 'analyte_name'] = fill_value
            else:
                analytes_index += 1
                fill_value = df['analyte_name'][analytes[analytes_index]].split(':')[1].strip()
        return df

    def process_data(self):
        current_dir = self.get_current_dir()
        file_type = self.file_type.get()
        files = get_files(current_dir, file_type)
        machine_type = self.machine_type.get()
        file_type = file_type[1:]
        ret_df_area = pd.DataFrame()
        ret_df_quantity = pd.DataFrame()
        if len(files):
            data_file = DataFile()
            for file in files:
                df = data_file.read(file, file_type)
                if machine_type == 'bruker':
                    df = janitor.clean_names(df, remove_special=True, case_type='snake')
                    df = df[BRUKER_VARIABLES]
                    df['quantity_units'] = pd.to_numeric(df['quantity_units'], errors='coerce')
                    df_area = df.pivot_table(index=['data_set', 'sample_type'], columns='analyte_name', values='area_of_pi')  # , aggfunc=np.mean)
                    df_quantity = df.pivot_table(index=['data_set', 'sample_type'], columns='analyte_name', values='quantity_units')  # , aggfunc=np.mean)
                else:
                    headers = df.loc[5].values.flatten().tolist()  # get the 5th row as headers
                    headers[0] = 'analyte_name'
                    headers[1] = 'hash'
                    df.columns = headers
                    df = janitor.clean_names(df, remove_special=True, case_type='snake')

                    headers = df.columns.to_list()
                    headers = [x.strip('_') for x in headers]  # remove leading and trailing '_' in variables
                    df.columns = headers

                    df.dropna(subset=['analyte_name'], inplace=True)  # drop empty rows 
                    df.reset_index(drop=True, inplace=True)  # reindex after dropping rows

                    df = self.fill_analyte_name(df)  # Compound: tryptophan occurs only once, fill it in rows below it
                    df = df[WATERS_VARIABLES]
                    df['conc'] = pd.to_numeric(df['conc'], errors='coerce')
                    
                    df["area"] = pd.to_numeric(df["area"], errors='coerce')
                    df_area = df.pivot_table(index=['sample_text', 'type'], columns='analyte_name', values='area')  # , fill_value=0)  # , aggfunc=np.mean)
                    df_quantity = df.pivot_table(index=['sample_text', 'type'], columns='analyte_name', values='conc')  # , fill_value=0)  # , aggfunc=np.mean)

                ret_df_area = ret_df_area.append(df_area)
                ret_df_quantity = ret_df_quantity.append(df_quantity)
        return ret_df_area, ret_df_quantity

    def long_to_wide(self):
        self.selected_cwd = True

        area_df, quantity_df = self.process_data()

        current_dir = self.get_current_dir()
        today = datetime.today()
        timestamp = f"{today.year}{today.month:02}{today.day:02}_{today.hour:02}{today.minute:02}{today.second:02}"
        filename = f"Flipped_{timestamp}.xlsx"
        output_path = os.path.join(current_dir, filename)

        machine_type = self.machine_type.get()
        sheet_name1 = 'Area'
        sheet_name2 = 'Conc'
        if machine_type == 'bruker':
            sheet_name1 = 'Area of PI'
            sheet_name2 = 'Quantity Units'

        with ExcelWriter(output_path) as writer:
            area_df.to_excel(writer, sheet_name1)
            quantity_df.to_excel(writer, sheet_name2)
            writer.save()

        self.process_message.configure(text=f"\nSaved as: {output_path}", fg="#006600", bg="#b3cccc")

    def select_cwd(self):
        old = self.config["cwd"]
        new = open_dir(old)
        self.selected_cwd = True
        if new:
            set_config("cwd", new)
            self.config["cwd"] = new
            self.cwd_label_text = new
            self.process_message.configure(text="\nFolder Selected")
            if old != new:
                self.cwd_label.configure(text="Changed Folder.", fg="green")
                self.dir_lf.configure(text=new, fg="#000")
            else:
                self.cwd_label.configure(text="", fg="blue")
                self.dir_lf.configure(text=new, fg="#000")
        elif old:
            self.dir_lf.configure(text=old, fg="#000")
            self.process_message.configure(text="\nSelected last Folder")
        self.process_button.configure(state=tk.NORMAL)

    def change_color(self):
        current_fg = self.dir_lf.cget("foreground")
        other = "#aaa"
        if not self.selected_cwd and current_fg in ["red", other]:
            next_fg = other if current_fg == "red" else "red"
            self.dir_lf.configure(fg=next_fg)
            root.after(600, self.change_color)

    def add_process_controls(self):
        process_lf = tk.LabelFrame(self.master, text="", padx=2, pady=2, relief=tk.FLAT, bg="#b3cccc")
        process_lf.pack(side=tk.BOTTOM, padx=2, pady=2)

        help_lf = tk.LabelFrame(process_lf, text="Background", padx=2, pady=2, relief=tk.SUNKEN, bg="#b3cccc")
        help_lf.pack(side=tk.TOP, padx=2, pady=(2, 30))

        help_message = tk.Label(help_lf, bg="#b3cccc", fg="midnightblue", font=("Arial", 10), justify=tk.LEFT)
        help_message.pack(side=tk.TOP, padx=2, pady=2)

        help_text = f"Please copy your files to a local folder on this PC and select that folder.\n\n"

        help_text += f"Columns in long format xlsx export files:\n\n"
        help_text += f"Bruker: {BRUKER_VARIABLES}\n\n"
        help_text += f"Waters: {WATERS_HELP_VARIABLES}\nFor Waters, analyte names appear as separate rows like Compound: tryptophan etc."
        help_message.configure(text=help_text)

        self.process_button = tk.Button(process_lf, text="Flip (Long to Wide)", command=self.long_to_wide, state=tk.DISABLED, font=("Arial", 12))
        self.process_button.pack(side=tk.TOP, padx=2, pady=2)

        self.process_message = tk.Label(process_lf, text=None, fg="green", bg="#b3cccc", font=("Arial", 10))
        self.process_message.pack(side=tk.BOTTOM, padx=2, pady=2, fill='both')

    def add_controls(self):
        cwd = self.config["cwd"]
        self.dir_lf = tk.LabelFrame(self.cwd_lf, text='Select your data directory', padx=2, pady=2, relief=tk.RIDGE, bg="#b3cccc", fg="red")
        self.dir_lf.pack(side=tk.LEFT, padx=8, pady=2)
        cwd_button = tk.Button(self.dir_lf, text="Select your Folder", command=self.select_cwd)
        cwd_button.pack(side=tk.LEFT, padx=2, pady=2)
        self.cwd_label = tk.Label(self.dir_lf, text=self.cwd_label_text, fg="#000", bg="#b3cccc")
        self.cwd_label.pack(side=tk.LEFT, padx=2, pady=2)

        type_lf = tk.LabelFrame(self.cwd_lf, text="File Type", padx=2, pady=2, relief=tk.RIDGE, bg="#b3cccc")
        type_lf.pack(side=tk.LEFT, padx=8, pady=2)
        types = [
            (".xlsx", ".xlsx"),
            (".csv", ".csv"),
            (".txt", ".txt")
        ]
        for text, ftype in types:
            tk.Radiobutton(type_lf, text=text, variable=self.file_type, value=ftype).pack(side=tk.LEFT, padx=2, pady=5)

        machine_lf = tk.LabelFrame(self.cwd_lf, text="Machine", padx=2, pady=2, relief=tk.RIDGE, bg="#b3cccc")
        machine_lf.pack(side=tk.LEFT, padx=8, pady=2)
        machines = [
            ("Bruker", "bruker"),
            ("Waters", "waters"),
        ]                
        for name, code in machines:
            tk.Radiobutton(machine_lf, text=name, variable=self.machine_type, value=code).pack(side=tk.LEFT, padx=2, pady=5)

        self.add_process_controls()

    def exit(self):
        clear_output()
        self.master.destroy()

    def add_options_bar(self):
        tk.Label(self, text="Flippy", bg="#b3cccc", font=("Arial", 16)).pack(side=tk.TOP, fill=tk.X)

        self.optionsbar = tk.Frame(self, bd=1, relief=tk.RIDGE, bg="#b3cccc")
        self.optionsbar.pack(side=tk.TOP, fill=tk.X)

        self.cwd_lf = tk.Frame(self.optionsbar, bg="#b3cccc", padx=4, pady=2, relief=tk.FLAT)
        self.cwd_lf.pack(side=tk.LEFT, padx=2, pady=2)

        self.add_controls()

    def create_widgets(self):
        self.menubar = tk.Menu(master=self.master, bg="#aaa")
        self.master.config(menu=self.menubar)
        self.menubar.add_command(label="Exit", command=self.exit, activebackground="#b3cccc", foreground="red")
        self.add_options_bar()
        self.process_message.configure(text=f"\nLast used directory: {self.config['cwd']}", fg="#666")

        self.change_color()


root = tk.Tk()
app = Application(root)
root.geometry("800x400")
root.configure(background='#b3cccc')
root.mainloop()
