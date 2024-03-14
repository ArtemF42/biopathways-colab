import pandas as pd
# from core import Dataset

db = pd.read_csv('./ligand_receptor/human_lr_pair.txt', sep='\t')


# def annotate(df: Dataset) -> Dataset:
#     genes = df.index.to_list()
#     pass


print(db)
