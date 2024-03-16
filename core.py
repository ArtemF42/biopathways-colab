from typing import List

import pandas as pd
import seaborn as sns


class PathwayDatabase:
    DATABASES = {
        'hallmark': 'h.all.v2023.2.Hs.symbols.gmt',
        'canonical pathways': 'c2.cp.v2023.2.Hs.symbols.gmt',
        'reactome': 'c2.cp.reactome.v2023.2.Hs.symbols.gmt',
        'transcription factor targets': 'c3.tft.v2023.2.Hs.symbols.gmt'
    }

    def __init__(self, database: str, pathway: str=None, genes: List[str]=None) -> None:
        self.pathways = {}
        
        if database != 'custom':
            with open(f'./biopathways/databases/{self.DATABASES[database]}') as file:
                for line in file:
                    pathway, _, *genes = line.removesuffix('\n').split('\t')
                    self.pathways[pathway] = genes
        else:
            self.pathways[pathway] = genes

    def __getitem__(self, pathway: str) -> List[str]:
        if pathway in self.pathways:
            return self.pathways[pathway]
    
    def search(self, query: str) -> List[str]:
        return [pathway for pathway in self.pathways if query in pathway]


def read(filepath: str) -> pd.DataFrame:
    match filepath.split('.')[-1]:
        case 'csv':
            sep = ','
        case 'tsv':
            sep = '\t'
        case _:
            raise Exception
        
    return pd.read_csv(filepath, sep=sep, index_col=0)


def significance(p_value: float) -> str:
    if p_value < 0.001:
        return '***'
    elif p_value < 0.01:
        return '**'
    elif p_value < 0.05:
        return '*'
    else:
        return ''


def heatmap(data: pd.DataFrame, annot: pd.DataFrame, method: str, metric: str):
    n_rows, n_cols = data.shape

    return sns.clustermap(data=data,
                          method=method,
                          metric=metric,
                          figsize=(n_cols * 2, n_rows * 0.5),
                          cmap='bwr',
                          center=0,
                          annot=annot.applymap(significance),
                          fmt='',
                          col_cluster=False)
