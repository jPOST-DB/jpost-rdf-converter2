
from __future__ import annotations

from dataclasses import dataclass, field

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from .protein import Protein
    from .peptide import Peptide
    from .dataset import DataSet

@dataclass
class Group:
    id: str | None = None
    proteins: list[Protein] = field(default_factory=list)
    peptides: list[Peptide] = field(default_factory=list)

    def __init__(self, id):
        self.id = id
        self.proteins = []
        self.peptides = []

    def get_id(self) -> str | None:
        return self.id

    def get_proteins(self) -> list[Protein]:
        return self.proteins
    
    def get_peptides(self) -> list[Peptide]:
        return self.peptides
    
    def add_protein(self, protein: Protein) -> None:
        found = False
        for p in self.proteins:
            if p.get_id() == protein.get_id():
                found = True

        if not found:
            self.proteins.append(protein)


    def add_peptide(self, peptide: Peptide) -> None:
        found = False
        for p in self.peptides:
            if p.get_dummy() == peptide.get_dummy():
                found = True

        if not found:
            self.peptides.append(peptide)

    def __str__(self) -> str:
        return f'Group(proteins={len(self.proteins)}, peptides={len(self.peptides)})'
    
    def __hash__(self) -> int:
        return id(self)
    
    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Group):
            return False
        return self is other    


    def to_ttl(self, f) -> None:
        f.write(f':{self.get_id()} a bid:ProteinGroup .\n\n')


    @staticmethod
    def create_groups(dataset: DataSet, proteins: list[Protein], optimized_proteins: list[Protein]) -> list[Group]:
        groups: list[Group] = []
        group_map: dict[str, Group] = {}

        for protein in optimized_proteins:
            uniprot = protein.get_uniprot()
            index = uniprot.rfind('-')
            if index >= 0:
                uniprot = uniprot[:index]

            group = protein.get_group()
            if group is None:
                group_id = f'PG{dataset.get_number()}_{len(groups)+1}'
                group = Group(group_id)
                groups.append(group)

            if uniprot not in group_map:
                group_map[uniprot] = group


        peptide_map: dict[Group, set[str]] = {}
        for protein in proteins:
            uniprot = protein.get_uniprot()
            index = uniprot.rfind('-')
            if index >= 0:
                uniprot = uniprot[:index]

            if uniprot in group_map:
                group = group_map[uniprot]
                protein.set_group(group)
                group.add_protein(protein)

                if group not in peptide_map:
                    peptide_map[group] = set()

                peptide_set = peptide_map[group]

                peptides = protein.search_peptides()
                for peptide in peptides:
                    peptide_set.add(peptide.get_dummy())
            else:
                protein.set_group(None)


        for protein in proteins:
            group = protein.get_group()
            if group is None:
                peptides = protein.search_peptides()
                
                for current in peptide_map.keys():
                    peptide_set = peptide_map.get(current)
                    flag = True
                    for peptide in peptides:
                        if peptide.get_dummy() not in peptide_set:
                            flag = False

                    if flag:
                        group = current

                if group is None:
                    id = f'PG{dataset.get_number()}_{len(groups)+1}'
                    group = Group(id)
                    peptide_map[group] = set()
                    for peptide in peptides:
                        peptide_map[group].add(peptide.get_dummy())
                    groups.append(group)

                protein.set_group(group)
        return groups
    
    @staticmethod
    def save_groups(f, groups: list[Group]) -> None:
        headers = ['Group ID', 'UniProt', 'Protein ID', 'Isoform', 'Protein Type', 'Leading Protein ID']
        f.write('\t'.join(headers) + '\n')

        for group in groups:
            for protein in group.get_proteins():
                row = [
                    group.get_id(), protein.get_uniprot(), protein.get_id()
                ]

                isoforms = ''
                for isoform in protein.get_isoforms():
                    if isoforms:
                        isoforms += ', '
                    isoforms += isoform.get_id()
                row.append(isoforms)

                types = []
                if protein.is_leading():
                    types.append('leading protein')
                if protein.is_anchor():
                    types.append('anchor protein')
                if protein.is_subset():
                    types.append('subset protein')
                if protein.is_same():
                    types.append('shared protein')
                row.append(', '.join(types))

                leading = ''
                for p in protein.get_leading_proteins():
                    if leading:
                        leading += ', '
                    leading += p.get_id()
                row.append(leading)

                f.write('\t'.join(row) + '\n')
