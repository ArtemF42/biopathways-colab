import pandas as pd

db = pd.read_csv('./ligand_receptor/human_lr_pair.txt', sep='\t')

print(db.head())
