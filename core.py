from typing import List

import pandas as pd
import seaborn as sns


from typing import List


class PathwayDatabase:
    DATABASES = {
        'hallmark': 'h.all.v2023.2.Hs.symbols.gmt',
        'canonical pathways': 'c2.cp.v2023.2.Hs.symbols.gmt',
        'reactome': 'c2.cp.reactome.v2023.2.Hs.symbols.gmt',
        'transcription factor targets': 'c3.tft.v2023.2.Hs.symbols.gmt'
    }

    def __init__(self, pathways: dict) -> None:
        self.pathways = pathways

    def __getitem__(self, pathway: str) -> List[str] | None:
        if pathway in self.pathways:
            return self.pathways[pathway]
    
    def search(self, query: str, by_gene: bool = False) -> List[str]:
        if by_gene:
            return [pathway for pathway in self.pathways if query in self.pathways[pathway]]
        else:
            return [pathway for pathway in self.pathways if query in pathway]
    
    @classmethod
    def load(cls, path: str, database: str):
        pathways = {}

        with open(f'{path}/{cls.DATABASES[database]}') as file:
            for line in file:
                pathway, _, *genes = line.removesuffix('\n').split('\t')
                pathways[pathway] = genes
        
        return cls(pathways)

    @classmethod
    @property
    def collection(cls) -> List:
        return list(cls.DATABASES)


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


def heatmap(data: pd.DataFrame, annot: pd.DataFrame | None, method: str, metric: str, vertical: bool):
    n_rows, n_cols = data.shape
    data.index.name = ''

    return sns.clustermap(data=data,
                          method=method,
                          metric=metric,
                          cbar_kws={'label': 'log2FoldChange'},
                          figsize=((n_cols * 2) if n_cols > 2 else (n_cols * 3), n_rows * 0.5),
                          dendrogram_ratio=(.1, .2),
                          cmap='bwr',
                          center=0 if annot is not None else None,
                          cbar_pos=(-0.1, .3, .04, .3),
                          annot=annot.applymap(significance) if annot is not None else None,
                          fmt='',
                          z_score=0 if annot is None else None,
                          col_cluster=vertical)
