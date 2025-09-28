from __future__ import annotations

from dataclasses import dataclass, field

from typing import TYPE_CHECKING
if TYPE_CHECKING:
     from .protein import Protein
     from .peptide import PeptideMatch


@dataclass
class Isoform:
    id: str | None = None
    uniprot: str | None = None
    type: str | None = None
    protein: Protein | None = None
    in_optimization_list: bool = False
    peptide_matchs: list[PeptideMatch] = field(default_factory=list)

    def __init__(self):
        self.id = None
        self.uniprot = None
        self.type = None
        self.protein = None
        self.in_optimization_list = False
        self.peptide_matchs = []

    def get_id(self) -> str | None:
        return self.id
    
    def set_id(self, id: str) -> None:
        self.id = id

    def get_uniprot(self) -> str | None:
        return self.uniprot
    
    def set_uniprot(self, uniprot: str) -> None:
        self.uniprot = uniprot

    def get_type(self) -> str | None:
        return self.type

    def set_type(self, type: str) -> None:
        self.type = type

    def get_protein(self) -> Protein | None:
        return self.protein
    
    def set_protein(self, protein: Protein) -> None:
        self.protein = protein

    def is_in_optimization_list(self) -> bool:
        return self.in_optimization_list

    def set_in_optimization_list(self, in_optimization_list: bool) -> None:
        self.in_optimization_list = in_optimization_list

    def get_peptide_matches(self) -> list[PeptideMatch]:
        return self.peptide_matchs
    
    def to_ttl(self, f) -> None:
        f.write(f':{self.protein.get_id()} jpost:hasIsoform :{self.id} .\n')
        f.write(f':{self.id}\n')
        f.write(f'    dct:identifier "{self.id}" ;\n')
        f.write(f'    rdfs:label "{self.uniprot}" ;\n')
        f.write(f'    rdfs:seeAlso isoforms:{self.uniprot} ;\n')
        f.write(f'    jpost:hasDatabaseSequence isoforms:{self.uniprot} ;\n')
        
        if self.is_in_optimization_list():
            f.write(f'    a jpost:RepresentativeIsoform ;\n')

            for peptide_match in self.peptide_matchs:
                peptide_match.to_ttl(f)

        f.write(f'    a jpost:ProteinIsoform .\n\n')
