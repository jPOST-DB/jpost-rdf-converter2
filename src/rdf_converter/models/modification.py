from __future__ import annotations

from dataclasses import dataclass, field
from ..utils.string_tool import is_not_empty
import xml.etree.ElementTree as ET
import requests 
import re
import os
from dotenv import load_dotenv
import logging

logger = logging.getLogger(__name__)

# unimod.xml の specificity classification -> jPOST ontology (jpost.owl JPO_022-039)
CLASSIFICATION_TO_JPO = {
    'Post-translational': 'JPO_022',
    'Co-translational': 'JPO_023',
    'Pre-translational': 'JPO_024',
    'Chemical derivative': 'JPO_025',
    'Artefact': 'JPO_026',
    'N-linked glycosylation': 'JPO_027',
    'O-linked glycosylation': 'JPO_028',
    'Other glycosylation': 'JPO_029',
    'Synth. pep. protect. gp.': 'JPO_030',
    'Isotopic label': 'JPO_031',
    'Non-standard residue': 'JPO_032',
    'Multiple': 'JPO_033',
    'Other': 'JPO_034',
    'AA substitution': 'JPO_035',
    'Cross-link': 'JPO_036',
    'CID cleavable cross-link': 'JPO_037',
    'Photo cleavable cross-link': 'JPO_038',
    'Other cleavable cross-link': 'JPO_039',
}

_unimod_specificity = None


def _load_unimod_specificity() -> dict:
    """UNIMOD_XML (未設定なら ./resources/unimod.xml) を読み、{record_id: [(site, position, classification), ...]} を返す"""
    global _unimod_specificity
    if _unimod_specificity is not None:
        return _unimod_specificity
    _unimod_specificity = {}
    load_dotenv()
    path = os.getenv('UNIMOD_XML', os.path.join('resources', 'unimod.xml'))
    if not os.path.exists(path):
        logger.warning(f'unimod.xml not found ({path}); jpost:modificationClass is not output')
        return _unimod_specificity
    ns = '{http://www.unimod.org/xmlns/schema/unimod_2}'
    for _, elem in ET.iterparse(path):
        if elem.tag == f'{ns}mod':
            specs = [(sp.get('site'), sp.get('position'), sp.get('classification'))
                     for sp in elem.findall(f'{ns}specificity')]
            _unimod_specificity[elem.get('record_id')] = specs
            elem.clear()
    return _unimod_specificity


def get_modification_class(unimod: str | None, site: str | None, title: str | None = None) -> str | None:
    """UNIMOD ID と修飾部位から modificationClass (JPO_xxx) を求める。
    site は 'C' / 'N-term' / 'Protein N-term' 等。title に '(Protein N-term)' があれば Protein N-term の分類を使う"""
    if unimod is None:
        return None
    specs = _load_unimod_specificity().get(str(unimod))
    if not specs:
        return None
    site = (site or '').strip()
    position = None
    m = re.match(r'(?i)^protein\s*([NC])-term$', site)
    if m:
        site, position = f'{m.group(1).upper()}-term', f'Protein {m.group(1).upper()}-term'
    elif site in ('N-term', 'C-term') and title and f'Protein {site}' in title:
        position = f'Protein {site}'

    def lookup(residue):
        candidates = [sp for sp in specs if sp[0] == residue]
        if position is not None:
            candidates = [sp for sp in candidates if sp[1] == position] or candidates
        if not candidates:
            return None
        # 部位が同じものが複数ある場合 (Any N-term / Protein N-term 等) は Any/Anywhere を優先
        candidates.sort(key=lambda sp: 0 if (sp[1] or '').startswith('Any') else 1)
        return CLASSIFICATION_TO_JPO.get(candidates[0][2])

    if site in ('N-term', 'C-term') or len(site) <= 1:
        clazz = lookup(site)
        if clazz is None and position is not None:
            # 例: UNIMOD:766 は site="M" position="Protein N-term" のみ定義
            classes = {CLASSIFICATION_TO_JPO.get(sp[2]) for sp in specs if sp[1] == position}
            clazz = classes.pop() if len(classes) == 1 else None
        return clazz
    # 複数残基表記 (例: NQ) は全残基で同じ分類のときのみ
    classes = {lookup(r) for r in site}
    return classes.pop() if len(classes) == 1 else None


@dataclass
class Modification:
    title: str | None = None
    unimod: str | None = None
    site: str | None = None
    clazz: str | None = None

    def __init__(self):
        self.title = None
        self.unimod = None
        self.site = None
        self.clazz = None

    def get_title(self) -> str | None:
        return self.title
    
    def set_title(self, title: str) -> None:
        self.title = title

    def get_unimod(self) -> str | None:
        return self.unimod
    
    def set_unimod(self, unimod: str) -> None:
        self.unimod = unimod

    def get_site(self) -> str | None:
        return self.site
    
    def set_site(self, site: str) -> None:
        self.site = site

    def get_class(self) -> str | None:
        return self.clazz
    
    def set_class(self, clazz: str) -> None:
        self.clazz = clazz

    def to_ttl(self, f) -> None:
        f.write(f'        a unimod:UNIMOD_{self.unimod} ;\n')
        if is_not_empty(self.site):
            f.write(f'        jpost:modificationSite "{self.site}" ;\n')
        clazz = self.clazz or get_modification_class(self.unimod, self.site, self.title)
        if clazz is not None:
            f.write(f'        jpost:modificationClass jpost:{clazz} ;\n')


    @staticmethod
    def get_modifications_from_jpost_repo(project_id: str):
        load_dotenv()
        repository_url = os.getenv('REPOSITORY_URL', 'https://repository.jpostdb.org/xml/')

        xml_url = f'{repository_url}{project_id}.0.xml'
        response = requests.get(xml_url)
        fixed_mods = []
        variable_mods = []

        pattern_fixed = re.compile(r"^fixed_mod\[(\d+)\]")
        pattern_variable = re.compile(r"^variable_mod\[(\d+)\]")
        pattern_fixed_site = re.compile(r"^fixed_mod\[(\d+)\]-site")
        pattern_variable_site = re.compile(r"^variable_mod\[(\d+)\]-site")

        if response.status_code == 200:
            root = ET.fromstring(response.content)

            file_tags = root.findall('FileList/File')
            result_file = None
            for file_tag in file_tags:
                type = file_tag.find('Type').text
                if type == 'result':
                    result_file = file_tag.find('Name').text

            if result_file is not None:
                result_url = f'https://storage.jpostdb.org/{project_id}/{result_file}'
                response = requests.get(result_url)
                fixed_mod_map = {}
                variable_mod_map = {}
                for line in response.iter_lines():
                    tokens = line.decode('utf-8').split('\t')
                    if len(tokens) >= 3:
                        if tokens[0] == 'MTD':
                            match_fixed = pattern_fixed.match(tokens[1])
                            match_variable = pattern_variable.match(tokens[1])
                            match_fixed_site = pattern_fixed_site.match(tokens[1])
                            match_variable_site = pattern_variable_site.match(tokens[1])

                            if match_fixed_site:
                                index = match_fixed_site.group(1)
                                site = tokens[2].strip()
                                if index in fixed_mod_map:
                                    fixed_mod_map[index]['site'] = site
                            elif match_variable_site:
                                index = match_variable_site.group(1)
                                site = tokens[2].strip()
                                if index in variable_mod_map:
                                    variable_mod_map[index]['site'] = site
                            elif (match_fixed or match_variable) and tokens[1].endswith('-position'):
                                # 'Protein N-term' 等は title に反映する (例: Acetyl (Protein N-term))
                                index = (match_fixed or match_variable).group(1)
                                mod_map = fixed_mod_map if match_fixed else variable_mod_map
                                if index in mod_map:
                                    mod_map[index]['position'] = tokens[2].strip()
                            elif match_fixed:
                                if not tokens[1].endswith('-position'):
                                    index = match_fixed.group(1)
                                    values = tokens[2].replace('[', '').replace(']', '').split(',')
                                    fixed_mod_map[index] = {'unimod': values[1].strip().replace('UNIMOD:', ''),
                                                            'name': values[2].strip()}
                            elif match_variable:
                                if not tokens[1].endswith('-position'):
                                    index = match_variable.group(1)
                                    values = tokens[2].replace('[', '').replace(']', '').split(',')
                                    variable_mod_map[index] = {'unimod': values[1].strip().replace('UNIMOD:', ''),
                                                            'name': values[2].strip()}
                            

                for mod_dic in fixed_mod_map.values():
                    modification = Modification()
                    modification.set_unimod(mod_dic['unimod'])
                    title = mod_dic['name']
                    if 'site' in mod_dic:
                        if mod_dic.get('position', '').startswith('Protein ') and mod_dic['site'] in ('N-term', 'C-term'):
                            title += f" ({mod_dic['position']})"
                        else:
                            title += f" ({mod_dic['site']})"
                        modification.set_site(mod_dic['site'])
                    modification.set_title(title)
                    fixed_mods.append(modification)

                for mod_dic in variable_mod_map.values():
                    modification = Modification()
                    modification.set_unimod(mod_dic['unimod'])
                    title = mod_dic['name']
                    if 'site' in mod_dic:
                        if mod_dic.get('position', '').startswith('Protein ') and mod_dic['site'] in ('N-term', 'C-term'):
                            title += f" ({mod_dic['position']})"
                        else:
                            title += f" ({mod_dic['site']})"
                        modification.set_site(mod_dic['site'])
                    modification.set_title(title)
                    variable_mods.append(modification)

        return fixed_mods, variable_mods
                
                        






