import os
import pandas as pd


class DataFile:
    @classmethod
    def read(cls, file, file_type):
        f = "read_" + file_type
        # for invalid file type, fn is the lambda function that returns error message
        fn = getattr(cls, f, lambda: "Invalid file type")
        return fn(file)

    def read_xlsx(file):
        return pd.read_excel(file)

    def read_TXT(file):
        with open(file) as tsv, open(file + ".csv", "w+") as tmp:
            first_row = True
            for aline in tsv:
                row = aline.strip('\n').split('\t')
                if first_row:
                    for i in  range(13 - len(row)):
                        row.append(f"Dummy{i}")  # a dummy header row
                    first_row = False
                else:
                    for i in  range(13 - len(row)):
                        row.append("")
                tmp.write(",".join(row) + "\n")
        df = pd.read_csv(file + ".csv")
        os.remove(file + ".csv")
        return df
