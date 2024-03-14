from typing import List


class PathwayDatabase:
    DATABASES = {
        'hallmark': 'h.all.v2023.2.Hs.symbols.gmt',
        'canonical pathways': 'c2.cp.v2023.2.Hs.symbols.gmt',
        'reactome': 'c2.cp.reactome.v2023.2.Hs.symbols.gmt',
        'transcription factor targets': 'c3.tft.v2023.2.Hs.symbols.gmt'
    }

    def __init__(self, database: str) -> None:
        self.pathways = {}

        with open(f'./biopathways/databases/{self.DATABASES[database]}') as file:
            for line in file:
                pathway, _, *genes = line.removesuffix('\n').split('\t')
                self.pathways[pathway] = genes

    def __getitem__(self, pathway: str) -> List[str]:
        if pathway in self.pathways:
            return self.pathways[pathway]
    
    def search(self, query: str) -> List[str]:
        return [pathway for pathway in self.pathways if query in pathway]
