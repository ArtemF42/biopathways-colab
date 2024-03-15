import pandas as pd
import numpy as np
from typing import Tuple

LRDB = pd.read_csv('./biopathways/ligand_receptor/human_lr_pair.txt', sep='\t')

LIGANDS = LRDB['ligand_gene_symbol'].unique()
RECEPTORS = LRDB['receptor_gene_symbol'].unique()


def annotate(df: pd.DataFrame, pval: pd.DataFrame, lr_choice: str) -> Tuple[pd.DataFrame, pd.DataFrame]:    
    match lr_choice:
        case 'receptors':
            df = df.loc[df.index.isin(RECEPTORS)]
        case 'ligands':
            df = df.loc[df.index.isin(LIGANDS)]
        case 'both':
            df = df.loc[np.concatenate([LIGANDS, RECEPTORS])]
            df.index = df.index.map(lambda gene: gene + ' (receptor)' 
                                    if gene in RECEPTORS else gene + ' (ligand)')
        case _:
            raise ValueError
    
    return df, pval.loc[df.index]
