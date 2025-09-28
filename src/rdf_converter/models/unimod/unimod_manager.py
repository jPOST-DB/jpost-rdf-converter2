from __future__ import annotations

from dataclasses import dataclass, field
import logging
from os import name
from dotenv import load_dotenv
import os
import xml.etree.ElementTree as ET


@dataclass
class UnimodElement:
    id: str
    site: str
    clazz: str
    title: str

    def __init__(self, id: str, site: str, clazz: str, title: str):
        self.id = id
        self.site = site
        self.clazz = clazz
        self.title = title

    def get_id(self) -> str:
        return self.id
    
    def set_id(self, id: str) -> None:  
        self.id = id
    
    def get_site(self) -> str:
        return self.site
    
    def set_site(self, site: str) -> None:  
        self.site = site
    
    def get_class(self) -> str:
        return self.clazz
    
    def set_class(self, clazz: str) -> None:    
        self.clazz = clazz
    
    def get_title(self) -> str:
        return self.title
    
    def set_title(self, title: str) -> None:    
        self.title = title
    
    def __str__(self):
        return f'UnimodElement(id={self.id}, site={self.site}, class={self.clazz}, title={self.title})' 


@dataclass
class TitleAndSite:
    title: str
    site: str

    def __init__(self, title: str, site: str):
        self.title = title
        self.site = site

    def get_title(self) -> str:
        return self.title
    
    def get_site(self) -> str:
        return self.site
    
    def __str__(self):
        return f'TitleAndSite(title={self.title}, site={self.site})'    
    

@dataclass
class UnimodManager:
    __instance = None

    @staticmethod
    def get_instance() -> UnimodManager:
        if UnimodManager.__instance is None:
            UnimodManager()
        return UnimodManager.__instance
    
    def __init__(self):
        if UnimodManager.__instance is not None:
            raise Exception("This class is a singleton!")
        else:            
            UnimodManager.__instance = self
            load_dotenv()
            unimod_xml = os.getenv('UNIMOD_XML', './unimod.xml')
            self.unimod_xml = unimod_xml
            self.list = UnimodManager.load_unimod_list(self.unimod_xml)



    def add(self, id, title, site, clazz) -> None:
        element = UnimodElement(id, site, clazz, title)
        self.list.append(element)


    def get_title_and_site(self, name: str) -> TitleAndSite:
        title = ''
        site = ''

        index1 = name.find('(')
        index2 = name.find(')')

        if index2 > index1 and index1 >= -1:
            title = name[:index1].strip()
            site = name[index1+1:index2].strip()
        
        return TitleAndSite(title, site)
    

    def search_element(self, id: str, site: str, is_protein: bool) -> UnimodElement | None:
        title = id.replace('UNIMOD_', '')
        index = title.find('#')
        if index > 0:
            title = title[:index]
        title = title.strip()

        if site.find('Protein') >= 0:
            site = site.replace('Protein', '').strip()
            is_protein = True

        if is_protein:
            site = f'Protein {site}'

        element = None
        for e in self.list:
            if e.get_id() == title:
                if site == '':
                    if element is None:
                        element = e
                elif e.get_site() == site:
                    element = e
        return element
    

    def search_element_by_name(self, name) -> UnimodElement | None:
        element = None

        title_and_site = UnimodManager.get_title_and_site(name)
        title = title_and_site.get_title()
        site = title_and_site.get_site()

        if site.find('Protein') >= 0:
            site = site.replace('Protein', '').strip()

        for e in self.list:
            if e.get_title() == title and e.get_site() == site:
                element = e

        return element


    def __str__(self):
        return f'UnimodManager(unimod_xml={self.unimod_xml}, list_size={len(self.list)})'
    

    @staticmethod
    def load_unimod_list(unimod_xml: str) -> list[UnimodElement]:
        list: list[UnimodElement] = []

        if not os.path.exists(unimod_xml):
            raise FileNotFoundError(f'Unimod XML file not found: {unimod_xml}')

        tree = ET.parse(unimod_xml)
        root = tree.getroot()

        ns = {"umod": "http://www.unimod.org/xmlns/schema/unimod_2"}

        mods = root.findall('umod:modifications/umod:mod', ns)

        for mod in mods:
            id = mod.get('record_id')
            title = mod.get('title')

            specificities = mod.findall('umod:specificity', ns)

            for specificity in specificities:
                clazz = specificity.get('classification')
                site = specificity.get('site')
                position = specificity.get('position')

                if position.lower().find('protein') >= 0 \
                        and (site.lower() == 'c-term' or site.lower() == 'n-term'):
                    site = f'Protein {site}'

                element = UnimodElement(id, site, clazz, title)
                list.append(element)

        return list
