import pandas as pd

# Write some information about the samples
df_full = pd.read_excel("mmc5.xlsx", header = 5, usecols="A:DE")
samples_info = df_full.iloc[5:, 0:4]
samples_info.to_csv('sample_information.csv', header=False, index = False)


filename = "structure_input.txt"

#%%
df = df_full.iloc[5:, 5:]
df.insert(loc=0, column = "ID", value = df_full.iloc[5:, 3].to_list())
df.insert(loc=1, column = "POP", value = df_full.iloc[5:, 1].to_list())
df.insert(loc=2, column = "LOC", value = df_full.iloc[5:, 2].to_list())
df.set_index("ID", inplace=True)
df.dropna(how='all', axis=1, inplace=True)

# Extract only European and Middle East samples
# df = df[[("EURO" in i or "MIDDLE" in i) for i in df["POP"]]]

pops = list(set(df["POP"]))
df["POP"] = [pops.index(i) for i in df["POP"]]
locs = list(set(df["LOC"]))
df["LOC"] = [locs.index(i) for i in df["LOC"]]


#%%

# Codierungsregel: A -> 0, C -> 1, G -> 2, T -> 3
code_map = {'A': '0', 'C': '1', 'G': '2', 'T': '3', 'N': '-1'}

# Write structure input file
def get_data(row):
    first_row = []
    second_row = []
    for value in row:
        value = str(value)
        if pd.notna(value) and isinstance(value, str):
            encoded_values = [code_map.get(base, 'N') for base in value if base in code_map]
            if len(encoded_values) == 2:
                first_row.append(encoded_values[0])
                second_row.append(encoded_values[1])
            else:
                first_row.append('-1')
                second_row.append('-1')
    return first_row, second_row

with open(filename, "w") as f:
    f.writelines('\t'.join(df.columns.to_list()[2:]) + '\n')
    for index, row in df.iterrows():
        data = row[2:].to_list()
        first_row, second_row = get_data(data)
        f.writelines('\t'.join([index, str(row["POP"]), str(row["LOC"])] + first_row) + '\n')
        f.writelines('\t'.join([index, str(row["POP"]), str(row["LOC"])] + second_row) + '\n')

print("The following information is necessary for the mainparms structure file:")        
print(df.shape)
