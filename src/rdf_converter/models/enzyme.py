from __future__ import annotations


from dataclasses import dataclass, field
import logging

import xml.etree.ElementTree as ET

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .dataset import DataSet

from .modification import Modification
from ..utils.string_tool import is_not_empty


logger = logging.getLogger(__name__)


@dataclass
class Enzyme:
    dataset: DataSet
    id: str | None = None
    enzymes: list[str] = field(default_factory=list)
    species: str | None = None
    fixed_mods: list[Modification] = field(default_factory=list)
    variable_mods: list[Modification] = field(default_factory=list)

    def __init__(self, dataset: DataSet):
        self.dataset = dataset
        dataset.set_enzyme(self)
        self.id = f'ENZ{dataset.get_number()}'
        self.species = None
        self.enzymes = []
        self.fixed_mods = []
        self.variable_mods = []


    def get_dataset(self) -> DataSet:
        return self.dataset
    
    def get_id(self) -> str | None:
        return self.id

    def get_species(self) -> str | None:
        return self.species
    
    def set_species(self, species: str) -> None:
        self.species = species

    def get_enzymes(self) -> list[str]:
        return self.enzymes
    
    def get_fixed_mods(self) -> list[Modification]:
        return self.fixed_mods
    
    def get_variable_mods(self) -> list[Modification]:
        return self.variable_mods
    
    def __str__(self):
        return f'Enzyme(id={self.id}, species={self.species}, enzymes={self.enzymes}, fixed_mods={self.fixed_mods}, variable_mods={self.variable_mods})'



    def to_ttl(self, f) -> None:
        enzyme = self
        dataset = self.get_dataset()

        f.write(f'bid:PRF{dataset.get_number()} jpost:hasEnzyme bid:{enzyme.get_id()} .\n')
        f.write(f'bid:{enzyme.get_id()}\n')

        for element in enzyme.get_enzymes():
            f.write(f'    jpost:enzyme obo:{element.replace(":", "_")} ;\n')


        for mod in enzyme.get_fixed_mods():
            f.write('    jpost:fixedModification [\n')
            mod.to_ttl(f)
            f.write('    ] ;\n')
            
        for mod in enzyme.get_variable_mods():
            f.write('    jpost:variableModification [\n')
            mod.to_ttl(f)
            f.write('    ] ;\n')

        f.write('    a jpost:EnzymeAndModifications .\n\n')


    @staticmethod
    def read_enzyme(dataset: DataSet, meta_path: str) -> Enzyme:
        tree = ET.parse(meta_path)
        root = tree.getroot()

        enzyme = Enzyme(dataset)

        tag = root.find('FileList/File/Profile/Enzyme_Mod')

        if tag is not None:
            taxonomy = tag.find('taxonomy')
            if taxonomy is not None and taxonomy.text is not None:
                enzyme.set_species(taxonomy.text.strip())
            
            enzyme_tag = tag.find('enzyme')
            if enzyme_tag is not None:
                enzyme_id = enzyme_tag.get('id')
                if is_not_empty(enzyme_id):
                    enzyme.get_enzymes().append(enzyme_id.strip())
            
            fixedMods, variableMods = Modification.get_modifications_from_jpost_repo(dataset.get_project().get_id())
            
            for mod in fixedMods:
                enzyme.fixed_mods.append(mod)

            for mod in variableMods:
                enzyme.variable_mods.append(mod)

        return enzyme