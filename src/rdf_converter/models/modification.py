from __future__ import annotations

from dataclasses import dataclass, field
from ..utils.string_tool import is_not_empty
import xml.etree.ElementTree as ET
import requests 
import re
import os
from dotenv import load_dotenv

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
        if is_not_empty(self.site):
            f.write(f'        jpost:modificationSite "{self.site}" ;\n')        
        f.write(f'        a unimod:UNIMOD_{self.unimod} \n')


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
                                index = match_fixed_site.group(0)
                                site = tokens[2].strip()
                                if index in fixed_mod_map:
                                    fixed_mod_map[index]['site'] = site
                            elif match_variable_site:
                                index = match_variable_site.group(0)
                                site = tokens[2].strip()
                                if index in variable_mod_map:
                                    variable_mod_map[index]['site'] = site
                            elif match_fixed:
                                if not tokens[1].endswith('-position'):
                                    index = match_fixed.group(0)
                                    values = tokens[2].replace('[', '').replace(']', '').split(',')
                                    fixed_mod_map[index] = {'unimod': values[1].strip().replace('UNIMOD:', ''),
                                                            'name': values[2].strip()}
                            elif match_variable:
                                if not tokens[1].endswith('-position'):
                                    index = match_variable.group(0)
                                    values = tokens[2].replace('[', '').replace(']', '').split(',')
                                    variable_mod_map[index] = {'unimod': values[1].strip().replace('UNIMOD:', ''),
                                                            'name': values[2].strip()}
                            

                for mod_dic in fixed_mod_map.values():
                    modification = Modification()
                    modification.set_unimod(mod_dic['unimod'])
                    title = mod_dic['name']
                    if 'site' in mod_dic:
                        title += f" ({mod_dic['site']})"
                        modification.set_site(mod_dic['site'])
                    modification.set_title(title)
                    fixed_mods.append(modification)

                for mod_dic in variable_mod_map.values():
                    modification = Modification()
                    modification.set_unimod(mod_dic['unimod'])
                    title = mod_dic['name']
                    if 'site' in mod_dic:
                        title += f" ({mod_dic['site']})"
                        modification.set_site(mod_dic['site'])
                    modification.set_title(title)
                    variable_mods.append(modification)

        return fixed_mods, variable_mods
                
                        






