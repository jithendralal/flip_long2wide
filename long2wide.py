# coding: utf-8

import pandas as pd
import janitor
import tkinter as tk
from datetime import datetime
from pandas import ExcelWriter
from mol_weights import *
from utils import *
from models import DataFile


class Application(tk.Frame):
    def __init__(self, master=None):
        tk.Frame.__init__(self, master)
        self.selected_cwd = False
        self.pack(fill='both', padx=2, pady=2)
        self.master.title('ANPC - Flippy')
        self.machine_type = tk.StringVar(value="Bruker")
        self.analysis_type = tk.StringVar(value="Amino Acids")
        self.load_config()
        self.create_widgets()
        self.df = pd.DataFrame()

    def load_config(self):
        self.config = get_config()

    def get_current_dir(self):
        return self.config["cwd"]

    def disable_analysis_types(self):
        for w in self.analysis_lf.winfo_children():
            w.configure(state=tk.DISABLED)

    def show_bruker_analysis_types(self):
        self.analysis_lf.configure(text="Bruker")
        self.disable_analysis_types()
        for w in self.analysis_lf.winfo_children():
            if w['text'] in self.bruker_analysis_types:
                w.configure(state=tk.NORMAL)

    def show_waters_analysis_types(self):
        self.analysis_lf.configure(text="Waters")
        self.disable_analysis_types()
        for w in self.analysis_lf.winfo_children():
            if w['text'] in self.waters_analysis_types:
                w.configure(state=tk.NORMAL)

    def show_sciex_analysis_types(self):
        self.analysis_lf.configure(text="Sciex")
        self.disable_analysis_types()
        for w in self.analysis_lf.winfo_children():
            if w['text'] in self.sciex_analysis_types:
                w.configure(state=tk.NORMAL)

    def toggle_unit_conc(self, enable):
        if enable:
            self.unit_conc_lf.configure(text="Tryptophan")
            for w in self.unit_conc_lf.winfo_children():
                w.configure(state=tk.NORMAL)
        else:
            self.unit_conc_lf.configure(text=" ")
            for w in self.unit_conc_lf.winfo_children():
                w.configure(state=tk.DISABLED)

    def toggle_double_conc(self, enable):
        if enable:
            self.double_conc_lf.configure(text="Amino Acids")
            for w in self.double_conc_lf.winfo_children():
                w.configure(state=tk.NORMAL)
        else:
            self.double_conc_lf.configure(text=" ")
            for w in self.double_conc_lf.winfo_children():
                w.configure(state=tk.DISABLED)

    def set_processing_options(self):
        analysis_type = self.get_analysis_type()
        if analysis_type == "Tryptophan":
            self.toggle_unit_conc(True)
            self.toggle_double_conc(False)
        elif analysis_type == "Amino Acids":
            self.toggle_double_conc(True)
            self.toggle_unit_conc(False)
        else:
            self.toggle_double_conc(False)
            self.toggle_unit_conc(False)

    def set_selections_text(self):
        machine = self.machine_type.get()
        if machine == 'Waters':
            self.show_waters_analysis_types()
        if machine == 'Sciex':
            self.show_sciex_analysis_types()
        if machine == 'Bruker':
            self.show_bruker_analysis_types()
        analysis_type = self.get_analysis_type()
        message = f"Folder: {self.get_current_dir()}\n"
        message += f"Machine: {machine}\n"
        message += f"File type: .{self.get_file_type()}\n"
        message += f"Analysis: {analysis_type}\n"
        self.selections_text.configure(text=message)
        self.set_processing_options()

    def fill_analyte_name(self, df):
        analytes = df[df['analyte_name'].apply(lambda x: str(x).startswith('Compound'))].index.tolist()
        analytes_index = 0
        fill_value = df['analyte_name'][analytes[analytes_index]].split(':')[1].strip()
        for i in range(2, len(df.index)):
            if not str(df['analyte_name'][i]).startswith('Compound'):
                df.at[i, 'analyte_name'] = fill_value
            else:
                analytes_index += 1
                fill_value = df['analyte_name'][analytes[analytes_index]].split(':')[1].strip()
        return df

    def extra_process_bruker(self, df_quantity):
        if self.double_conc.get():
            for col in df_quantity.columns.to_list():
                df_quantity[col] = df_quantity[col].apply(double_it)
        return df_quantity

    def process_bruker(self, df):
        df = janitor.clean_names(df, remove_special=True, case_type='snake')
        df = df[BRUKER_VARIABLES]
        df['quantity_units'] = pd.to_numeric(df['quantity_units'], errors='coerce')
        df_area = df.pivot_table(index=['data_set', 'sample_type'], columns='analyte_name', values='area_of_pi')  # , aggfunc=np.mean)
        df_quantity = df.pivot_table(index=['data_set', 'sample_type'], columns='analyte_name', values='quantity_units')  # , aggfunc=np.mean)
        df_rt = df.pivot_table(index=['data_set', 'sample_type'], columns='analyte_name', values='rt_min')  # , aggfunc=np.mean)
        df_quantity = self.extra_process_bruker(df_quantity)
        return df_area, df_quantity, df_rt

    def extra_process_waters(self, df_quantity):
        if self.unit_conc.get():
            for col in df_quantity.columns.to_list():
                mol_weight = mol_weights[col]
                df_quantity[col] = df_quantity[col].apply(unit_conc, args=(mol_weight,))
        return df_quantity

    def process_waters(self, df):
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
        df["rt"] = pd.to_numeric(df["rt"], errors='coerce')
        df_area = df.pivot_table(index=['sample_text', 'type'], columns='analyte_name', values='area')  # , fill_value=0)  # , aggfunc=np.mean)
        df_quantity = df.pivot_table(index=['sample_text', 'type'], columns='analyte_name', values='conc')  # , fill_value=0)  # , aggfunc=np.mean)
        df_rt = df.pivot_table(index=['sample_text', 'type'], columns='analyte_name', values='rt')  # , fill_value=0)  # , aggfunc=np.mean)

        df_quantity = self.extra_process_waters(df_quantity)

        return df_area, df_quantity, df_rt

    def get_file_type(self):
        ret = 'TXT'
        if self.machine_type.get() == 'Bruker':
            ret = 'xlsx'
        return ret
    
    def get_analysis_type(self):
        return self.analysis_type.get()

    def machine_change(self):
        self.analysis_type.set("test")
        self.set_selections_text()

    def process_files(self):
        current_dir = self.get_current_dir()
        machine_type = self.machine_type.get()
        file_type = self.get_file_type()
        files = get_files(current_dir, "."+file_type)

        sheet_name1 = 'Area'
        sheet_name2 = 'Conc'
        sheet_name3 = 'RT'
        if machine_type == 'Bruker':
            sheet_name1 = 'Area of PI'
            sheet_name2 = 'Quantity Units'

        today = datetime.today()
        timestamp = f"{today.year}{today.month:02}{today.day:02}_{today.hour:02}{today.minute:02}{today.second:02}"

        if len(files):
            data_file = DataFile()
            for file in files:
                if '_flipped' not in file:
                    df = data_file.read(file, file_type)
                    if machine_type == 'Bruker':
                        df_area, df_quantity, df_rt = self.process_bruker(df)
                    else:
                        df_area, df_quantity, df_rt = self.process_waters(df)

                    out_filename = f"{file}_{timestamp}_flipped.xlsx"
                    output_path = os.path.join(current_dir, out_filename)
                    with ExcelWriter(output_path) as writer:
                        df_area.to_excel(writer, sheet_name1)
                        df_quantity.to_excel(writer, sheet_name2)
                        df_rt.to_excel(writer, sheet_name3)
                        writer.save()

    def long_to_wide(self):
        self.selected_cwd = True
        self.process_files()
        self.process_message.configure(text=f"Completed.", fg="#006600", bg="#ddd")

    def select_cwd(self):
        old = self.config["cwd"]
        new = open_dir(old)
        self.selected_cwd = True
        if new:
            set_config("cwd", new)
            self.config["cwd"] = new
            self.process_message.configure(text="New folder selected")
            if old != new:
                self.dir_lf.configure(text="", fg="#000")
            else:
                self.dir_lf.configure(text="", fg="#000")
        elif old:
            self.dir_lf.configure(text="", fg="#000")
            self.process_message.configure(text=f"Selected same folder: {old}")
        for w in self.machine_lf.winfo_children():
            w.configure(state=tk.NORMAL)
        self.process_button.configure(state=tk.NORMAL)
        self.set_selections_text()

    def change_color(self):
        current_fg = self.dir_lf.cget("foreground")
        other = "#aaa"
        if not self.selected_cwd and current_fg in ["red", other]:
            next_fg = other if current_fg == "red" else "red"
            self.dir_lf.configure(fg=next_fg)
            root.after(600, self.change_color)

    def add_feedback_controls(self):
        feedback_lf = tk.LabelFrame(self.master, text="Your selections:", fg='#666', padx=2, pady=2, relief=tk.FLAT, bg="#b3cccc")
        feedback_lf.pack(side=tk.TOP, padx=2, pady=2)

        self.selections_text = tk.Label(feedback_lf, bg="#b3cccc", fg="midnightblue", font=("Arial", 11), justify=tk.LEFT)
        self.selections_text.pack(side=tk.TOP, padx=2, pady=2)

    def add_process_controls(self):
        process_lf = tk.LabelFrame(self.master, text="", padx=2, pady=2, relief=tk.FLAT, bg="#b3cccc")
        process_lf.pack(side=tk.BOTTOM, padx=0, pady=0, fill=tk.X)

        self.process_button = tk.Button(process_lf, text="Flip", activebackground='palegreen', width=20,
                                        command=self.long_to_wide, state=tk.DISABLED, font=("Arial", 12))
        self.process_button.pack(side=tk.TOP, padx=2, pady=(0, 20))

        self.process_message = tk.Label(process_lf, fg="green", bg="#ddd", font=("Arial", 10))
        self.process_message.pack(side=tk.BOTTOM, padx=0, pady=0, fill=tk.X)

    def add_help(self):
        help_lf = tk.LabelFrame(self.cwd_lf, text="", padx=2, pady=2, relief=tk.FLAT, bg="#ccc")
        help_lf.pack(side=tk.TOP, padx=1, pady=(1,10))

        help_text = f"Please copy your files to an empty folder and select that folder (Each file will be processed separately)\n\n"
        help_text += f"Files should have the columns:\n"
        help_text += f"Bruker (xlsx): {BRUKER_VARIABLES}\n"
        help_text += f"Waters (TXT) : {WATERS_HELP_VARIABLES} (analyte names appear on separate lines, eg. Compound: tryptophan){' '*32}"

        help_message = tk.Label(help_lf, bg="#ccc", fg="midnightblue", font=("Arial", 10), justify=tk.LEFT, text=help_text)
        help_message.pack(side=tk.TOP, padx=10, pady=0)

    def _add_machine_selector(self):
        self.machine_lf = tk.LabelFrame(self.cwd_lf, text="Machine", padx=2, pady=2, relief=tk.FLAT, bg="#ccc")
        self.machine_lf.pack(side=tk.LEFT, padx=8, pady=2)
        machines = [
            ("Bruker", "Bruker"),
            ("Waters", "Waters"),
            ("Sciex", "Sciex"),
        ]                
        for name, code in machines:
            tk.Radiobutton(self.machine_lf, text=name, variable=self.machine_type, bd=0, command=self.machine_change,
                           activebackground='palegreen', state=tk.DISABLED,
                           value=code, relief=tk.SOLID).pack(anchor=tk.W, padx=2, pady=2)

    def _add_analysis_type(self):
        self.analysis_lf = tk.LabelFrame(self.cwd_lf, text=" ", padx=2, pady=2, relief=tk.FLAT, bg="#ccc")
        self.analysis_lf.pack(side=tk.LEFT, padx=8, pady=2)

        self.bruker_analysis_types = [
            "Amino Acids",
        ]
        self.sciex_analysis_types = [
            "Targeted Lipids",
            "Lipid Mediators",
        ]
        self.waters_analysis_types = [
            "Tryptophan",
            "Bile Acids",
            "Paracetamol",
            "SCFAs",
        ]
        analysis_types = self.bruker_analysis_types + self.waters_analysis_types + self.sciex_analysis_types
        for name in analysis_types:
            wat_rb = tk.Radiobutton(self.analysis_lf, text=name, variable=self.analysis_type, bd=0, command=self.set_selections_text,
                           activebackground='palegreen', state=tk.DISABLED,
                           value=name, relief=tk.SOLID).pack(anchor=tk.W, padx=2, pady=2)

    def _add_double_conc_selector(self):
        self.double_conc_lf = tk.LabelFrame(self.cwd_lf, text=" ", padx=2, pady=2, relief=tk.FLAT, bg="#ccc")
        self.double_conc_lf.pack(side=tk.LEFT, padx=8, pady=2)
        self.double_conc = tk.IntVar(value=1)
        tk.Checkbutton(self.double_conc_lf, text="Double Qty.", state=tk.DISABLED, variable=self.double_conc, bg="#ddd", activebackground="palegreen").pack(side=tk.LEFT, padx=2, pady=2)

    def _add_unit_conc_selector(self):
        self.unit_conc_lf = tk.LabelFrame(self.cwd_lf, text=" ", padx=2, pady=2, relief=tk.FLAT, bg="#ccc")
        self.unit_conc_lf.pack(side=tk.LEFT, padx=8, pady=2)
        self.unit_conc = tk.IntVar(value=1)
        tk.Checkbutton(self.unit_conc_lf, text="Conc. Unit (ng/mL to uM)", state=tk.DISABLED, variable=self.unit_conc, bg="#ddd", activebackground="palegreen").pack(side=tk.LEFT, padx=2, pady=2)

    def add_controls(self):
        cwd = self.config["cwd"]
        self.dir_lf = tk.LabelFrame(self.cwd_lf, text='Important !', padx=2, pady=2, relief=tk.FLAT, bg="#ccc", fg="red")
        self.dir_lf.pack(side=tk.LEFT, padx=8, pady=2)
        cwd_button = tk.Button(self.dir_lf, text="Select folder", command=self.select_cwd, activebackground='palegreen')
        cwd_button.pack(side=tk.LEFT, padx=2, pady=2)
        self._add_machine_selector()
        self._add_analysis_type()
        self._add_double_conc_selector()
        self._add_unit_conc_selector()
        self.add_process_controls()
        self.add_feedback_controls()

    def exit(self):
        self.master.destroy()
    
    def close(self, event):
        self.exit()

    def add_options_bar(self):
        tk.Label(self, text="FlipPy", bg="#b3cccc", font=("Arial", 16)).pack(side=tk.TOP, fill=tk.X)
        tk.Label(self, text="(Cast Long to Wide)", bg="#b3cccc", font=("Arial", 11)).pack(side=tk.TOP, fill=tk.X)

        self.optionsbar = tk.Frame(self, bd=1, relief=tk.FLAT, bg="#ddd")
        self.optionsbar.pack(side=tk.TOP, fill=tk.X)

        self.cwd_lf = tk.Frame(self.optionsbar, bg="#ddd", padx=4, pady=2, relief=tk.FLAT)
        self.cwd_lf.pack(side=tk.LEFT, fill=tk.X)
        
        self.add_help()
        self.add_controls()

    def create_widgets(self):
        self.menubar = tk.Menu(master=self.master, bg="#aaa")
        self.master.config(menu=self.menubar)
        self.menubar.add_command(label="Exit <Esc>", command=self.exit, activebackground="#b3cccc", foreground="red")
        self.add_options_bar()
        self.process_message.configure(text=f"Last used directory: {self.config['cwd']}", fg="#666")

        self.change_color()
        self.master.bind('<Escape>', self.close)


root = tk.Tk()
app = Application(root)
root.geometry("900x600")
root.configure(background='#b3cccc')
root.mainloop()
