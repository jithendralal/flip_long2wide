# coding: utf-8
import os
import pandas as pd
import janitor
from datetime import datetime
from pandas import ExcelWriter
from webbrowser import open as url_open

from utils import *
from models import DataFile


def get_file_type(analysis_type):
    file_type = 'TXT'
    if analysis_type == 'Bruker Amino Acids':
        file_type = 'xlsx'
    return file_type


def get_uploads_dir():
    return os.path.join(os.getcwd(), "uploads")


def load_compounds():
    if os.path.isfile(COMPOUNDS_FILE):
        with open(COMPOUNDS_FILE) as f:
            compounds = json.load(f)
    return compounds


def fill_analyte_name(df):
    analytes = df.loc[df['analyte_name'].str.startswith('Compound'), ['analyte_name']]['analyte_name'].tolist()
    analytes_index = 0
    fill_value = analytes[analytes_index].split(':')[1].strip()
    for i in range(2, len(df.index)):
        if not str(df['analyte_name'][i]).startswith('Compound'):
            df.at[i, 'analyte_name'] = fill_value
        else:
            analytes_index += 1
            fill_value = analytes[analytes_index].split(':')[1].strip()
    return df


def double_quantity_bruker(df_quantity):
    if double_conc:
        for col in df_quantity.columns.to_list():
            df_quantity[col] = df_quantity[col].apply(double_it)
    return df_quantity


def process_bruker(df):
    if analysis_type == 'Bruker Amino Acids':
        df = janitor.clean_names(df, remove_special=True, case_type='snake')
        try:
            df = df[BRUKER_VARIABLES]
        except Exception as e:
            print("wrong parameters")
            return 'wrong parameters', str(e)
        df['quantity_units'] = pd.to_numeric(df['quantity_units'], errors='coerce')
        df_area = df.pivot_table(index=['data_set', 'sample_type'], columns='analyte_name', values='area_of_pi')
        df_quantity = df.pivot_table(index=['data_set', 'sample_type'], columns='analyte_name',
                                     values='quantity_units')
        df_rt = df.pivot_table(index=['data_set', 'sample_type'], columns='analyte_name', values='rt_min')
        df_quantity = double_quantity_bruker(df_quantity)
        return df_area, df_quantity, df_rt


def unit_conversion_waters(df_quantity):
    try:
        if unit_conc:
            for col in df_quantity.columns.to_list():
                if col != "Sample ID":
                    mol_weight = compounds[col.lower()]
                    df_quantity[col] = df_quantity[col].apply(unit_conc, args=(mol_weight,))
        if unit_conc_micro:
            for col in df_quantity.columns.to_list():
                if col != "Sample ID":
                    mol_weight = compounds[col.lower()]
                    df_quantity[col] = df_quantity[col].apply(unit_conc_micro, args=(mol_weight,))
        return df_quantity, None
    except KeyError as ke:
        missing_compound = str(ke).strip("'")
        return "missing compound", missing_compound


def process_conversion(df_quantity):
    result = unit_conversion_waters(df_quantity)
    if isinstance(result[0], str):
        return "missing compound", result[1]
    return result


def process_waters(df):
    if analysis_type == 'Waters Tryptophan' or analysis_type == 'Waters Bile Acids':
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

            df = fill_analyte_name(df)  # Compound: tryptophan occurs only once, fill it in rows below it

            df = df[WATERS_VARIABLES]  # ---------Note: for any new columns update this in utils and use it below--
            df['conc'] = pd.to_numeric(df['conc'], errors='coerce')
            df["area"] = pd.to_numeric(df["area"], errors='coerce')
            df["rt"] = pd.to_numeric(df["rt"], errors='coerce')

            # ---------Note: for any new columns update WATERS_VARIABLES and use it below------
            df_area = df.pivot_table(index=['sample_text', 'type'], columns='analyte_name', values='area')
            df_quantity = df.pivot_table(index=['sample_text', 'type'], columns='analyte_name', values='conc')
            df_rt = df.pivot_table(index=['sample_text', 'type'], columns='analyte_name', values='rt')

            result = unit_conversion_waters(df_quantity)
            if isinstance(result[0], str):
                return "missing compound", result[1]

            df_grouped_conc = group_conc_waters(df)

            return df_area, df_quantity, df_rt, df_grouped_conc
        except Exception as e:
            return "wrong parameters", str(e)


def group_conc_waters(df):
    df1 = df.loc[df["sample_text"] != "Double Blank"]
    df1 = df1.sort_values(by=['analyte_name', 'sample_text'], ignore_index=True)
    an = None
    st = None
    suffix = 0
    for index, row in df1.iterrows():
        if row['analyte_name'] != an:
            suffix = 0
            an = row['analyte_name']
            st = None
        if row['sample_text'] != st:
            st = row['sample_text']
            suffix = 0
        else:
            suffix += 1

        if suffix > 0:
            df1.at[index, 'sample_text'] = row['sample_text'] + str(suffix)
    df_out = df1.pivot_table(index=['sample_text'], columns='analyte_name', values='conc')
    return df_out


def process_files():
    file_type = get_file_type(analysis_type)
    files = get_files(uploads_dir, "." + file_type)

    sheet_name1 = 'Area'
    sheet_name2 = 'Conc'
    sheet_name3 = 'RT'
    sheet_name4 = 'Conc Grouped'
    if analysis_type == 'Bruker Amino Acids':
        sheet_name1 = 'Area of PI'
        sheet_name2 = 'Quantity Units'

    today = datetime.today()
    timestamp = f"{today.year}{today.month:02}{today.day:02}_{today.hour:02}{today.minute:02}{today.second:02}"

    if len(files):
        data_file = DataFile()
        for file in files:
            if '_flipped' not in file:
                df_grouped = None
                df = data_file.read(file, file_type)
                if analysis_type == 'Bruker Amino Acids':
                    result = process_bruker(df)
                    if isinstance(result[0], str) and result[0] == 'wrong parameters':
                        return result
                    df_area, df_quantity, df_rt = result[0], result[1], result[2]
                elif analysis_type == 'Waters Conversion':
                    result = process_conversion(df)
                    if isinstance(result[0], str) and (result[0] == 'missing compound' or
                                                       result[0] == 'wrong parameters'):
                        return result
                    out_filename = f"{file}_{timestamp}_converted.xlsx"
                    output_path = os.path.join(uploads_dir, out_filename)
                    with ExcelWriter(output_path) as writer:
                        result[0].to_excel(writer, "converted")
                        writer.save()
                    return 'completed'
                else:
                    result = process_waters(df)
                    if isinstance(result[0], str) and (result[0] == 'missing compound' or
                                                       result[0] == 'wrong parameters'):
                        return result
                    df_area, df_quantity, df_rt, df_grouped = result[0], result[1], result[2], result[3]

                out_filename = f"{file}_{timestamp}_flipped.xlsx"
                output_path = os.path.join(uploads_dir, out_filename)
                with ExcelWriter(output_path) as writer:
                    df_area.to_excel(writer, sheet_name1)
                    df_quantity.to_excel(writer, sheet_name2)
                    df_rt.to_excel(writer, sheet_name3)
                    if len(result) == 4:
                        df_grouped.to_excel(writer, sheet_name4)
                    writer.save()
    else:
        return 'not found'
    return 'completed'


def long_to_wide():
    result = process_files()
    if result == 'completed':
        print(f"Completed processing.")
        print(f"\n\nPlease check the uploads folder")
    elif result == 'not found':
        print(f"Files of type .{get_file_type()} not found in the uploads directory")
        print(f"Please select correct options")
    elif isinstance(result[0], str) and result[0] == 'wrong parameters':
        print(f"Please check your selections / file structure.")
        print(f"\nLook for missing columns, if you have modified the exported file.")
    elif isinstance(result[0], str) and result[0] == 'missing compound':
        print(f"Missing compound: {result[1]}")
        get_mol_mass(result[1])


def validate_mass():
    if missing_compound_mass and missing_compound_mass.isnumeric() and float(missing_compound_mass) > 0:
        add_compound()
        long_to_wide()
    else:
        missing_compound_mass = ""


def add_compound():
    compounds[missing_compound] = float(missing_compound_mass)
    with open(COMPOUNDS_FILE, "w+") as f:
        json.dump(compounds, f, indent=4)
    return


def get_search_url(missing_compound):
    return SEARCH_URL + missing_compound


def get_mol_mass(missing_compound):
    print(missing_compound)
    url_open(get_search_url(missing_compound.replace('_', ' ')))


SEARCH_URL = "https://google.com.au/search?q=molar mass of "

IMPLEMENTED = [
    'Bruker Amino Acids',
    'Waters Tryptophan',
    'Waters Bile Acids',
    'Waters Conversion'
]

COMPOUNDS_FILE = "compounds.json"

missing_compound = ""
missing_compound_mass = ""

analysis_type = IMPLEMENTED[0]
df = pd.DataFrame()
compounds = load_compounds()

uploads_dir = get_uploads_dir()

double_conc = False
unit_conc = False
unit_conc_micro = False


def main():
    long_to_wide()


if __name__=="__main__":
    main()
