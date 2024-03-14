import pandas as pd

db = pd.read_csv('./ligand_receptor/human_lr_pair.txt', sep='\t')


def annotate(df: pd.DataFrame, lr_choice: str) -> pd.DataFrame:
    genes = df.index
    ligands = [genes in db['ligand_gene_symbol'].unique()]
    receptors = [genes in db['receptor_gene_symbol'].unique()]
    
    match lr_choice:
        case 'receptors':
            df = df.loc[ligands]
        case 'ligands':
            df = df.loc[receptors]
        case 'both':
            df = df.loc[set(ligands, receptors)]
            df.index = df.index.map(
                lambda gene: gene + ' (receptor)' if gene in receptors else gene + ' (ligand)'
            )
        case _:
            raise ValueError
    
    return df