from __future__ import annotations

from dataclasses import dataclass, field
from ..utils.string_tool import is_not_empty
import logging
import xml.etree.ElementTree as ET

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .dataset import DataSet



logger = logging.getLogger(__name__)

@dataclass
class SubFractionation:
    id: str
    fractionation: str
    fractions: str
    replicates: str
    type: str

    def get_id(self) -> str:
        return self.id
    
    def get_fractionation(self) -> str:
        return self.fractionation
    
    def get_fractions(self) -> str:
        return self.fractions
    
    def get_replicates(self) -> str:
        return self.replicates
    
    def get_type(self) -> str:
        return self.type


@dataclass
class Fractionation:
    dataset: DataSet
    id: str | None = None
    subcellular_fractionation: str | None = None
    subcellular_fractions: str | None = None
    subcellular_replicates: str | None = None
    peptide_fractionation: str | None = None
    peptide_fractions: str | None = None
    peptide_replicates: str | None = None
    protein_fractionation: str | None = None
    protein_fractions: str | None = None
    protein_replicates: str | None = None
    subfractionations: list[SubFractionation] | None = None


    def __init__(self, dataset: DataSet):
        self.dataset = dataset
        self.id = f'FRC{dataset.get_number()}'
        dataset.set_fractionation(self)


    def get_dataset(self) -> DataSet:
        return self.dataset
    
    def get_id(self) -> str | None:
        return self.id
    
    def get_subcellular_fractionation(self) -> str | None:
        return self.subcellular_fractionation
    
    def set_subcellular_fractionation(self, subcellular_fractionation: str) -> None:
        self.subcellular_fractionation = subcellular_fractionation

    def get_subcellular_fractions(self) -> str | None:
        return self.subcellular_fractions
    
    def set_subcellular_fractions(self, subcellular_fractions: str) -> None:
        self.subcellular_fractions = subcellular_fractions

    def get_subcellular_replicates(self) -> str | None:
        return self.subcellular_replicates

    def set_subcellular_replicates(self, subcellular_replicates: str) -> None:
        self.subcellular_replicates = subcellular_replicates

    def get_peptide_fractionation(self) -> str | None:
        return self.peptide_fractionation
    
    def set_peptide_fractionation(self, peptide_fractionation: str) -> None:
        self.peptide_fractionation = peptide_fractionation

    def get_peptide_fractions(self) -> str | None:
        return self.peptide_fractions
    
    def set_peptide_fractions(self, peptide_fractions: str) -> None:
        self.peptide_fractions = peptide_fractions

    def get_peptide_replicates(self) -> str | None:
        return self.peptide_replicates
    
    def set_peptide_replicates(self, peptide_replicates: str) -> None:
        self.peptide_replicates = peptide_replicates

    def get_protein_fractionation(self) -> str | None:
        return self.protein_fractionation
    
    def set_protein_fractionation(self, protein_fractionation: str) -> None:
        self.protein_fractionation = protein_fractionation

    def get_protein_fractions(self) -> str | None:
        return self.protein_fractions
    
    def set_protein_fractions(self, protein_fractions: str) -> None:
        self.protein_fractions = protein_fractions

    def get_protein_replicates(self) -> str | None:
        return self.protein_replicates
    
    def set_protein_replicates(self, protein_replicates: str) -> None:
        self.protein_replicates = protein_replicates

    def get_subfractionations(self) -> list[SubFractionation]:
        if self.subfractionations is None:
            self.subfractionations = self.create_subfractionations()
        return self.subfractionations

    def create_subfractionations(self) -> list[SubFractionation]:
        subfractionations = []

        counter = 0

        if is_not_empty(self.subcellular_fractionation) \
                or is_not_empty(self.subcellular_fractions) \
                or is_not_empty(self.subcellular_replicates):
            counter
            sub_fractionation = SubFractionation(
                f'{self.id}_{counter}',
                self.subcellular_fractionation,
                self.subcellular_fractions,
                self.subcellular_replicates,
                'subcellular'
            )
            subfractionations.append(sub_fractionation)

        if is_not_empty(self.peptide_fractionation) \
                or is_not_empty(self.peptide_fractions) \
                or is_not_empty(self.peptide_replicates):
            counter += 1
            sub_fractionation = SubFractionation(
                f'{self.id}_{counter}',
                self.peptide_fractionation,
                self.peptide_fractions,
                self.peptide_replicates,
                'peptide'
            )
            subfractionations.append(sub_fractionation) 

        if is_not_empty(self.protein_fractionation) \
                or is_not_empty(self.protein_fractions) \
                or is_not_empty(self.protein_replicates):
            counter += 1
            sub_fractionation = SubFractionation(
                f'{self.id}_{counter}',
                self.protein_fractionation,
                self.protein_fractions,
                self.protein_replicates,
                'protein'
            )
            subfractionations.append(sub_fractionation)

        return subfractionations

    def __str__(self) -> str:
        return f'Fractionation(id={self.id}, subcellular_fractionation={self.subcellular_fractionation}, subcellular_fractions={self.subcellular_fractions}, subcellular_replicates={self.subcellular_replicates}, peptide_fractionation={self.peptide_fractionation}, peptide_fractions={self.peptide_fractions}, peptide_replicates={self.peptide_replicates}, protein_fractionation={self.protein_fractionation}, protein_fractions={self.protein_fractions}, protein_replicates={self.protein_replicates})'


    def to_ttl(self, f) -> None:
        dataset = self.get_dataset()
        fractionation = self

        for sub_fractionation in fractionation.get_subfractionations():
            f.write(f'bid:PRF{dataset.get_number()} jpost:hasFractionation bid:{sub_fractionation.get_id()} .\n')
            f.write(f'bid:{sub_fractionation.get_id()} jpost:hasFractionType [\n')

            label = sub_fractionation.get_fractionation()
            if is_not_empty(label):
                f.write(f'        rdfs:label "{label}" ;\n')

            replicates = sub_fractionation.get_replicates()
            if is_not_empty(replicates):
                try:
                    int(replicates)
                    f.write(f'        jpost:replicates {replicates} ;\n')
                except ValueError:
                    pass

            fractions = sub_fractionation.get_fractions()
            if is_not_empty(fractions):
                try:
                    int(fractions)
                    f.write(f'        jpost:fractions {fractions} ;\n')
                except ValueError:
                    pass
            
            type = sub_fractionation.get_type()
            clazz = 'jpost:Fractionation'
            if type == 'subcellular':
                clazz = 'jpost:SubcellularFractionation'
            elif type == 'protein':
                clazz = 'jpost:ProteinFractionation'
            elif type == 'peptide':
                clazz = 'jpost:PeptideFractionation'

            f.write(f'        a {clazz}\n')
            f.write('    ] .\n')
            f.write('\n')


    @staticmethod
    def read_fractionation(dataset: DataSet, xml_file: str) -> Fractionation:
        tree = ET.parse(xml_file)
        root = tree.getroot()

        fractionation = Fractionation(dataset)

        fractionation_tag = root.find('FileList/File/Profile/Fractionation')
        if fractionation_tag is not None:
            subcellular = fractionation_tag.findtext('subcellular')
            if subcellular is not None:
                fractionation.set_subcellular_fractionation(subcellular)

            peptide = fractionation_tag.find('peptide')
            if peptide is not None:
                fractionation.set_peptide_fractionation(peptide.text)
                fractions = peptide.get('fraction')
                if fractions is not None:
                    fractionation.set_peptide_fractions(fractions)
                replicates = peptide.get('replicate')
                if replicates is not None:
                    fractionation.set_peptide_replicates(replicates)

            protein = fractionation_tag.find('protein')
            if protein is not None:
                fractionation.set_protein_fractionation(protein.text)
                fractions = protein.get('fraction')
                if fractions is not None:
                    fractionation.set_protein_fractions(fractions)
                protein_replicates = protein.get('replicate')
                if protein_replicates is not None:
                    fractionation.set_protein_replicates(protein_replicates)

        return fractionation

