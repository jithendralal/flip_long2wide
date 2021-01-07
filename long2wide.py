# coding: utf-8

import pandas as pd
import janitor
import tkinter as tk
from datetime import datetime
from pandas import ExcelWriter
from tkinter import scrolledtext
from mol_weights import *
from utils import *
from models import DataFile


IMPLEMENTED = [
    ('Bruker', 'Amino Acids'),
    ('Waters', 'Tryptophan'),
]


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
        self.config_window = False
        self.messagebox = None

    def load_config(self):
        self.config = get_config()

    def get_current_dir(self):
        return self.config["cwd"]

    def disable_analysis_types(self):
        for w in self.analysis_lf.winfo_children():
            w.configure(state=tk.DISABLED)

    def show_bruker_analysis_types(self):
        self.analysis_lf.configure(text="Analysis (Bruker)")
        self.disable_analysis_types()
        for w in self.analysis_lf.winfo_children():
            if w['text'] in self.bruker_analysis_types:
                w.configure(state=tk.NORMAL)

    def show_waters_analysis_types(self):
        self.analysis_lf.configure(text="Analysis (Waters)")
        self.disable_analysis_types()
        for w in self.analysis_lf.winfo_children():
            if w['text'] in self.waters_analysis_types:
                w.configure(state=tk.NORMAL)

    def show_sciex_analysis_types(self):
        self.analysis_lf.configure(text="Analysis (Sciex)")
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
        self.status_message.configure(text=f" ", fg="#000")
        machine = self.machine_type.get()
        if machine == 'Waters':
            self.show_waters_analysis_types()
        if machine == 'Sciex':
            self.show_sciex_analysis_types()
        if machine == 'Bruker':
            self.show_bruker_analysis_types()
        analysis_type = self.get_analysis_type()
        message = f"Folder: {self.get_current_dir()}"
        self.selections_folder_text.configure(text=message)
        message = f"Machine: {machine}"
        self.selections_machine_text.configure(text=message)
        message = f"File type: .{self.get_file_type()}"
        self.selections_file_text.configure(text=message)
        message = f"Analysis: {analysis_type}"
        self.selections_analysis_text.configure(text=message)
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

    def not_implemented(self):
        pass

    def process_bruker(self, df):
        analysis_type = self.analysis_type.get()
        if analysis_type == 'Amino Acids':
            df = janitor.clean_names(df, remove_special=True, case_type='snake')
            try:
                df = df[BRUKER_VARIABLES]
            except Exception:
                return 'wrong parameters'
            df['quantity_units'] = pd.to_numeric(df['quantity_units'], errors='coerce')
            df_area = df.pivot_table(index=['data_set', 'sample_type'], columns='analyte_name', values='area_of_pi')  # , aggfunc=np.mean)
            df_quantity = df.pivot_table(index=['data_set', 'sample_type'], columns='analyte_name', values='quantity_units')  # , aggfunc=np.mean)
            df_rt = df.pivot_table(index=['data_set', 'sample_type'], columns='analyte_name', values='rt_min')  # , aggfunc=np.mean)
            df_quantity = self.extra_process_bruker(df_quantity)
            return df_area, df_quantity, df_rt

    def extra_process_waters(self, df_quantity):
        if self.unit_conc.get():
            for col in df_quantity.columns.to_list():
                mol_weight = mol_weights[col.lower()]
                df_quantity[col] = df_quantity[col].apply(unit_conc, args=(mol_weight,))
        return df_quantity

    def process_waters(self, df):
        analysis_type = self.analysis_type.get()
        if analysis_type =='Tryptophan':
            try:
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
            except Exception:
                return 'wrong parameters'

    def get_file_type(self):
        file_type = 'TXT'
        if self.machine_type.get() == 'Bruker':
            file_type = 'xlsx'
        return file_type

    def get_analysis_type(self):
        return self.analysis_type.get()

    def machine_change(self):
        self.analysis_type.set(" ")
        self.set_selections_text()

    def process_files(self):
        current_dir = self.get_current_dir()
        machine_type = self.machine_type.get()
        analysis_type = self.analysis_type.get()
        if analysis_type == ' ':
            message = f"Analysis type not selected."
            self.status_message.configure(text=message, fg="red")
            self.show_messagebox("Error", message)
            return
        if (machine_type, analysis_type) not in IMPLEMENTED:
            message = f"{machine_type}-{analysis_type} not implemented."
            self.status_message.configure(text=message, fg="red")
            self.show_messagebox("Error", message)
            return

        file_type = self.get_file_type()
        files = get_files(current_dir, "." + file_type)

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
                        result = self.process_bruker(df)
                        if result == 'wrong parameters':
                            return result
                        df_area, df_quantity, df_rt = result[0], result[1], result[2]
                    else:
                        result = self.process_waters(df)
                        if result == 'wrong parameters':
                            return result
                        df_area, df_quantity, df_rt = result[0], result[1], result[2]

                    out_filename = f"{file}_{timestamp}_flipped.xlsx"
                    output_path = os.path.join(current_dir, out_filename)
                    with ExcelWriter(output_path) as writer:
                        df_area.to_excel(writer, sheet_name1)
                        df_quantity.to_excel(writer, sheet_name2)
                        df_rt.to_excel(writer, sheet_name3)
                        writer.save()
        else:
            return 'not found'
        return 'completed'

    def show_messagebox(self, title, message):
        if not self.messagebox:
            self.messagebox = tk.Toplevel(self)
            self.messagebox.wm_title(title)
            self.messagebox.title(title)
            self.messagebox_label = tk.Label(self.messagebox, text=message)
            self.messagebox_label.pack(padx=5, pady=5)
            self.messagebox_button = tk.Button(self.messagebox, text="OK", command=lambda: self.messagebox.withdraw())
            self.messagebox_button.pack(padx=5, pady=10)
            self.messagebox.wm_protocol("WM_DELETE_WINDOW", lambda: self.on_delete_child(self.messagebox))
        else:
            self.messagebox.attributes('-topmost', 'true')
            self.messagebox.deiconify()
            self.messagebox.wm_title(title)
            self.messagebox.title(title)
            self.messagebox_label.configure(text=message)

    def long_to_wide(self):
        self.status_message.configure(text=f"Processing...", fg="#000066", bg="#ddd")
        self.selected_cwd = True
        result = self.process_files()
        if result == 'completed':
            message1 = f"Completed processing."
            message2 = f"\n Please check the folder: {self.config['cwd']}"
            self.status_message.configure(text=message1, fg="#006600", bg="#ddd")
            self.show_messagebox("Success", message1 + message2)
        elif result == 'not found':
            message1 = f"Files of type .{self.get_file_type()} not found in the directory"
            message2 = f"\n{self.config['cwd']}\nPlease select correct options"
            self.status_message.configure(text=message1, fg="#ff0000", bg="#ddd")
            self.show_messagebox("Error", message1 + message2)
        elif result == 'wrong parameters':
            message1 = f"Please check your selections / file structure."
            message2 = f"\nLook for missing columns, if you have modified the exported file."
            self.status_message.configure(text=message1, fg="#ff0000", bg="#ddd")
            self.show_messagebox("Error", message1 + message2)

    def select_cwd(self):
        old = self.config["cwd"]
        new = open_dir(old)
        self.selected_cwd = True
        if new:
            set_config("cwd", new)
            self.config["cwd"] = new
            self.status_message.configure(text="New folder selected")
            if old != new:
                self.dir_lf.configure(text="", fg="#000")
            else:
                self.dir_lf.configure(text="", fg="#000")
        elif old:
            self.dir_lf.configure(text="", fg="#000")
            self.status_message.configure(text=f"Selected same folder: {old}")
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
        feedback_lf = tk.LabelFrame(self.master, text="Your selections:", fg='#444', bg="#ccc", relief=tk.FLAT)
        feedback_lf.pack(side=tk.TOP, pady=10)
        
        f = tk.Frame(feedback_lf, bg="#ddd")
        f.pack(pady=3)
        self.selections_folder_text = tk.Label(f, bg="#ddd", fg="#222", font=("Arial", 10), anchor="w")
        self.selections_folder_text.pack(side=tk.TOP, padx=2, pady=1, fill=tk.BOTH, expand=1)

        h = 15
        w = 180

        f = tk.Frame(feedback_lf, height=h, width=w, bg="#eee")
        f.pack()
        f.pack_propagate(0)
        self.selections_machine_text = tk.Label(f, bg="#eee", fg="#222", font=("Arial", 10), anchor="w")
        self.selections_machine_text.pack(side=tk.TOP, padx=2, pady=1, fill=tk.BOTH, expand=1)

        f = tk.Frame(feedback_lf, height=h, width=w, bg="#ddd")
        f.pack()
        f.pack_propagate(0)
        self.selections_file_text = tk.Label(f, bg="#ddd", fg="#222", font=("Arial", 10), anchor="w")
        self.selections_file_text.pack(side=tk.TOP, padx=2, pady=1, fill=tk.BOTH, expand=1)

        f = tk.Frame(feedback_lf, height=h, width=w, bg="#eee")
        f.pack()
        f.pack_propagate(0)
        self.selections_analysis_text = tk.Label(f, bg="#eee", fg="#222", font=("Arial", 10), anchor="w")
        self.selections_analysis_text.pack(side=tk.TOP, padx=2, pady=1, fill=tk.BOTH, expand=1)

    def add_process_controls(self):
        status_bar = tk.LabelFrame(self.master, text="", padx=2, pady=2, relief=tk.FLAT, bg="#b3cccc")
        status_bar.pack(side=tk.BOTTOM, padx=0, pady=0, fill=tk.X)

        self.process_button = tk.Button(status_bar, text="Flip", activebackground='palegreen', width=20,
                                        command=self.long_to_wide, state=tk.DISABLED, font=("Arial", 12))
        self.process_button.pack(side=tk.TOP, padx=2, pady=(0, 20))

        self.status_message = tk.Label(status_bar, fg="green", bg="#ddd", font=("Arial", 10))
        self.status_message.pack(side=tk.BOTTOM, padx=0, pady=0, fill=tk.X)

    def add_help(self):
        help_lf = tk.LabelFrame(self.cwd_lf, text="", padx=2, pady=2, relief=tk.FLAT, bg="#ccc")
        help_lf.pack(side=tk.TOP, padx=1, pady=(1,10))

        help_text = f"Please copy your files to an empty folder and select that folder (Each file will be processed separately)\n\n"
        help_text += f"Files should have the columns:\n"
        help_text += f"Bruker (xlsx): {BRUKER_VARIABLES}\n"
        help_text += f"Waters (TXT) : {WATERS_HELP_VARIABLES} (analyte names appear on separate lines, eg. Compound: tryptophan){' '*32}"

        help_message = tk.Label(help_lf, bg="#ccc", fg="#222", font=("Arial", 10), justify=tk.LEFT, text=help_text)
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
        self.analysis_lf = tk.LabelFrame(self.cwd_lf, text="Analysis", padx=2, pady=2, relief=tk.FLAT, bg="#ccc")
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
        self.dir_lf = tk.LabelFrame(self.cwd_lf, text='Select', padx=2, pady=2, relief=tk.FLAT, bg="#ccc", fg="red")
        self.dir_lf.pack(side=tk.LEFT, padx=8, pady=2)
        cwd_button = tk.Button(self.dir_lf, text="Folder", command=self.select_cwd, activebackground='palegreen')
        cwd_button.pack(side=tk.LEFT, padx=2, pady=2)
        self._add_machine_selector()
        self._add_analysis_type()
        self._add_double_conc_selector()
        self._add_unit_conc_selector()
        self.add_process_controls()
        self.add_feedback_controls()

    def exit(self):
        self.master.destroy()

    def on_delete_child(self, w):
        w.destroy()
        self.config_window = None

    def help(self, event):
        self.show_help()

    def show_help(self):
        if not self.config_window:
            self.config_window = tk.Toplevel(self)
            self.config_window.wm_title("Help")
            self.config_window.wm_protocol("WM_DELETE_WINDOW", lambda: self.on_delete_child(self.config_window))

            help_text = f"FlipPy Instructions:\n"
            tk.Label(self.config_window, text=help_text, bg="#ddd", font=("Arial", 12), anchor="w").pack(fill=tk.X)

            help_text = f"\nPlease copy your files to an empty folder and select that folder."
            help_text += f"\nEach file will be processed separately and each will be saved as an xlsx.\n"
            tk.Label(self.config_window, text=help_text, bg="#ddd", font=("Arial", 10), justify="left").pack(fill=tk.X)

            help_text = f"\nFiles should have these columns:\n\n"
            help_text += f"Bruker (xlsx): {BRUKER_VARIABLES}\n"
            help_text += f"Waters (TXT) : {WATERS_HELP_VARIABLES} (analyte names appear on separate lines, eg. Compound: tryptophan)\n"
            tk.Label(self.config_window, text=help_text, bg="#ddd", font=("Arial", 10), justify="left").pack()

            tk.Label(self.config_window, text="Molecular weights:", bg="#ddd", font=("Arial", 12), anchor="w").pack(fill=tk.X)

            ctext = scrolledtext.ScrolledText(self.config_window, height=20, font=('Courier', 10))
            ctext.pack(padx=3, pady=3, fill=tk.BOTH)
            mol_list = [f"{k: <50}- {mol_weights[k]}\n" for k in mol_weights.keys()]
            help_text = "".join(mol_list)
            ctext.insert("end", help_text)
        else:
            tk.raise_above_all(self.config_window)

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

        filemenu = tk.Menu(self.menubar, tearoff=0)
        filemenu.add_command(label=f"Help           F1", command=self.show_help, activebackground="palegreen")
        filemenu.add_separator()
        filemenu.add_command(label=f"Exit            Esc", command=self.exit, activebackground="palegreen")
        self.menubar.add_cascade(label="File", menu=filemenu)

        self.add_options_bar()
        self.status_message.configure(text=f"Last used folder: {self.config['cwd']}", fg="#666")

        self.change_color()
        self.master.bind('<Escape>', self.close)
        self.master.bind('<F1>', self.help)
        self.select_cwd()


root = tk.Tk()
app = Application(root)
root.geometry("900x600")
root.configure(background='#b3cccc')
root.mainloop()
