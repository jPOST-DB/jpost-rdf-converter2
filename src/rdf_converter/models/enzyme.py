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

        f.write(f'bid:PRF{dataset.get_number()} jpost:hasEnzymeAndModification bid:{enzyme.get_id()} .\n')
        f.write(f'bid:{enzyme.get_id()}\n')

        for element in enzyme.get_enzymes():
            f.write(f'    jpost:enzyme obo:{element.replace(":", "_")} ;\n')

        if is_not_empty(enzyme.get_species()):
            f.write(f'    jpost:species taxonomy:{enzyme.get_species()} ;\n')


        for mod in enzyme.get_fixed_mods():
            f.write('    jpost:fixedModification [\n')
            mod.to_ttl(f)
            f.write('    ] ;\n')
            
        for mod in enzyme.get_variable_mods():
            f.write('    jpost:variableModification [\n')
            mod.to_ttl(f)
            f.write('    ] ;\n')

        f.write('    a jpost:EnzymeAndModification .\n\n')


    @staticmethod
    def read_modifications(tag, name: str) -> list[Modification]:
        """<fixedModification id="UNIMOD:4#C">Carbamidomethyl (C)</fixedModification> 形式を読む"""
        mods = []
        titles = set()
        for mod_tag in tag.findall(name):
            mod_id = (mod_tag.get('id') or '').strip()
            title = (mod_tag.text or '').strip()
            if not mod_id.startswith('UNIMOD:') or not title or title in titles:
                continue
            titles.add(title)
            unimod, _, site = mod_id[len('UNIMOD:'):].partition('#')
            if site in ('ProteinN-term', 'Protein N-term'):
                site = 'N-term'
            elif site in ('ProteinC-term', 'Protein C-term'):
                site = 'C-term'
            modification = Modification()
            modification.set_unimod(unimod.strip())
            modification.set_title(title)
            if site:
                modification.set_site(site.strip())
            mods.append(modification)
        return mods

    @staticmethod
    def read_enzyme(dataset: DataSet, meta_path: str) -> Enzyme:
        tree = ET.parse(meta_path)
        root = tree.getroot()

        enzyme = Enzyme(dataset)

        tag = root.find('FileList/File/Profile/Enzyme_Mod')

        if tag is not None:
            taxonomy = tag.find('taxonomy')
            if taxonomy is not None and is_not_empty(taxonomy.get('id')):
                enzyme.set_species(taxonomy.get('id').strip())
            
            # <enzyme> は複数ありうる (id のないもの (例: chymotrypsin) もある)。id のあるものを全て出力する
            for enzyme_tag in tag.findall('enzyme'):
                enzyme_id = enzyme_tag.get('id')
                if is_not_empty(enzyme_id) and enzyme_id.strip() not in enzyme.get_enzymes():
                    enzyme.get_enzymes().append(enzyme_id.strip())
            
            fixedMods, variableMods = Modification.get_modifications_from_jpost_repo(dataset.get_project().get_id())
            if not fixedMods and not variableMods:
                # repository の result ファイルが mzTab でない (MGF 等) 場合は、メタデータ XML の修飾定義 (UNIMOD 付き) を使う
                fixedMods = Enzyme.read_modifications(tag, 'fixedModification')
                variableMods = Enzyme.read_modifications(tag, 'variableModification')
                if fixedMods or variableMods:
                    logger.info('Modifications are taken from the metadata XML (not found in repository mzTab)')
            
            for mod in fixedMods:
                enzyme.fixed_mods.append(mod)

            for mod in variableMods:
                enzyme.variable_mods.append(mod)

        return enzyme