from typing import List


class PathwayDatabase:
    def __init__(self, database: str) -> None:
        self.pathways = {}

        with open(f'./biopathways/databases/{database}.gmt') as file:
            for line in file:
                pathway, _, *genes = line.removesuffix('\n').split('\t')
                self.pathways[pathway] = genes

    def __getitem__(self, pathway: str) -> List[str]:
        if pathway in self.pathways:
            return self.pathways[pathway]
    
    def search(self, query: str) -> List[str]:
        return [pathway for pathway in self.pathways if query in pathway]
