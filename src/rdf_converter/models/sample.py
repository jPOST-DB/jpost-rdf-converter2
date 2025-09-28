from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional, List

import xml.etree.ElementTree as ET
import logging

from ..utils.string_tool import is_not_empty

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .dataset import DataSet

logger = logging.getLogger(__name__)


@dataclass
class Sample:
    dataset: DataSet
    id: str | None = None
    species: str | None = None
    sample_type: str | None = None
    cell_line_id: str | None = None
    cell_line: str | None = None
    disease: str | None = None
    disease_id: str | None = None
    disease_class: str | None = None
    organ: str | None = None

    def __init__(self, dataset: DataSet):
        self.dataset = dataset
        dataset.set_sample(self)
        self.id = f'SMP{dataset.get_number()}'
        self.species = None
        self.sample_type = None
        self.cell_line_id = None
        self.cell_line = None
        self.disease = None
        self.disease_id = None
        self.disease_class = None
        self.organ = None


    def get_dataset(self) -> DataSet:
        return self.dataset
    
    def get_id(self) -> str | None:
        return self.id
    
    def get_species(self) -> str | None:
        return self.species
    
    def set_species(self, species: str) -> None:
        self.species = species

    def get_sample_type(self) -> str | None:
        return self.sample_type
    
    def set_sample_type(self, sample_type: str) -> None:
        self.sample_type = sample_type

    def get_cell_line_id(self) -> str | None:
        return self.cell_line_id
    
    def set_cell_line_id(self, cell_line_id: str) -> None:
        self.cell_line_id = cell_line_id

    def get_cell_line(self) -> str | None:
        return self.cell_line
    
    def set_cell_line(self, cell_line: str) -> None:
        self.cell_line = cell_line

    def get_disease(self) -> str | None:
        return self.disease
    
    def set_disease(self, disease: str) -> None:
        self.disease = disease

    def get_disease_id(self) -> str | None:
        return self.disease_id
    
    def set_disease_id(self, disease_id: str) -> None:
        self.disease_id = disease_id

    def get_disease_class(self) -> str | None:
        return self.disease_class
    
    def set_disease_class(self, disease_class: str) -> None:
        self.disease_class = disease_class

    def get_organ(self) -> str | None:
        return self.organ
    
    def set_organ(self, organ: str) -> None:
        self.organ = organ

    def __str__(self) -> str:
        return f'Sample(id={self.id}, species={self.species}, sample_type={self.sample_type}, cell_line_id={self.cell_line_id}, cell_line={self.cell_line}, disease={self.disease}, disease_id={self.disease_id}, disease_class={self.disease_class}, organ={self.organ})'

    def to_ttl(self, f) -> None:
        sample = self
        dataset = self.get_dataset()
        
        f.write(f'bid:PRF{dataset.get_number()} jpost:hasSample bid:{sample.get_id()} .\n')
        f.write(f'bid:{sample.get_id()}\n')

        species = sample.get_species()
        sample_type = sample.get_sample_type()
        cell_line_id = sample.get_cell_line_id()
        cell_line = sample.get_cell_line()
        organ = sample.get_organ()
        disease_class = sample.get_disease_class()
        disease_id = sample.get_disease_id()
        disease = sample.get_disease()

        if is_not_empty(species):
            f.write(f'    jpost:species taxonomy:{species} ;\n')

        if is_not_empty(sample_type):
            f.write(f'    jpost:sampleType ncit:{sample_type} ;\n')


        if is_not_empty(cell_line_id):
            prefix = 'obo'
            if cell_line_id.startswith('JPO'):
                prefix = 'jpost'
            f.write(f'    jpost:cellLine {prefix}:{cell_line_id.replace(':', '_')} ;\n')
        elif is_not_empty(cell_line):
            f.write('    jpost:cellLine [\n')
            f.write(f'        rdfs:label "{cell_line}" ;\n')
            f.write('    ] ;\n')

        if is_not_empty(organ):
            f.write(f'    jpost:organ ncit:{organ} ;\n')

        if is_not_empty(disease_class):
            prefix = 'obo'
            if disease_class.startswith('JPO'):
                prefix = 'jpost'
            f.write(f'    jpost:diseaseClass {prefix}:{disease_class.replace(':', '_')} ;\n')

        if is_not_empty(disease_id):
            prefix = 'obo'
            if disease_id.startswith('JPO'):
                prefix = 'jpost'
            f.write(f'    jpost:disease {prefix}:{disease_id.replace(":", "_")} ;\n')
        elif is_not_empty(disease):
            f.write('    jpost:disease [\n')
            f.write(f'        rdfs:label "{disease}" ;\n')
            f.write('    ] ;\n')
        f.write('    a jpost:Sample .\n\n')        

    @staticmethod
    def read_sample(dataset: DataSet, xml_file: str) -> Sample:        
        tree = ET.parse(xml_file)
        root = tree.getroot()

        sample = Sample(dataset)

        species = root.find('presetSummary/Species/PresetElement')
        if species is not None:
            id = species.get('id')
            sample.set_species(id)  

        note = root.find('FileList/File/Profile/Sample/note')
        if note is not None and note.text is not None:
            lines = note.text.split('\n')
            for line in lines:
                tokens = line.split('|')
                if len(tokens) >= 3:
                    key = tokens[0].strip()
                    if key == 'Sample Type':
                        sample.set_sample_type(tokens[2].strip())
                    elif key == 'Cell line':
                        cell_line = tokens[1].strip()
                        cell_line_id = tokens[2].strip()
                        sample.set_cell_line(cell_line)
                        sample.set_cell_line_id(cell_line_id)
                    elif key == 'Disease class':
                        sample.set_disease_class(tokens[2].strip())
                    elif key == 'Disease':
                        sample.set_disease(tokens[2].strip())
                    elif key == 'Organ':
                        sample.set_organ(tokens[2].strip())

        return sample

