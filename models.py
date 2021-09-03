import os
import pandas as pd


class DataFile:
    @classmethod
    def read(cls, file, file_type):
        f = "read_" + file_type
        fn = getattr(cls, f, lambda: "Invalid file type")
        return fn(file)

    def read_xlsx(file):
        return pd.read_excel(file)

    def read_TXT(file):
        with open(file) as tsv, open(file + ".csv", "w+") as tmp:
            top_rows = []
            i = 0
            col_count = 0

            for aline in tsv:
                row = aline.strip('\n').replace(',',';').split('\t')  # if there is a comma in data replace it with semi-colon
                if i <= 6:
                    top_rows.append(row)

                if i == 6:  # this is the header row (7th)
                    col_count = len(row)  # no of actual columns in TXT
                    
                    for col in range(col_count-1):
                        top_rows[0].append(f"Dummy{col}")  # a dummy header row
                    
                    for j in range(1, 6):  # make the number of columns same for first 6 rows
                        top_row = top_rows[j]
                        for col in range(col_count - len(top_row)):
                            top_row.append("")

                    for row in top_rows:
                        tmp.write(",".join(row) + "\n")  # write to csv

                if i > 6:
                    for col in  range(col_count - len(row)):
                        row.append("")
                    tmp.write(",".join(row) + "\n")

                i = i + 1
        
        df = pd.read_csv(file + ".csv")
        os.remove(file + ".csv")
        return df
