from __future__ import annotations

from dataclasses import dataclass, field

from typing import Optional


from xml.etree import ElementTree as ET

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from .dataset import DataSet


@dataclass
class RawData:
    file: str
    id: str

    def __init__(self, file: str, id: str):
        self.file = file
        self.id = id
    
    def get_file(self) -> str:
        return self.file

    def get_id(self) -> str:
        return self.id


@dataclass
class RawDataList:
    dataset: DataSet
    files: list[RawData] = field(default_factory=list)
    counter: int = 0

    def __init__(self, dataset: DataSet):
        self.dataset = dataset
        dataset.set_rawdata_list(self)
        self.files = []
        self.counter = 0        

    
    def get_dataset(self) -> DataSet:
        return self.dataset
    
    def get_list(self) -> list[RawData]:
        return self.files
    
    def add_file(self, file: str):
        self.counter += 1
        id = self.dataset.get_id().replace('DS', 'RAW') + f'_{self.counter}'
        
        rawdata = RawData(file, id)
        self.files.append(rawdata)


    def get_file_id(self, file: str) -> Optional[str]:  
        id = None
        for rawdata in self.files:
            if rawdata.get_file() == file:
                id = rawdata.get_id()
        return id
        
    def __str__(self):
        return f'RawDataList(dataset={self.dataset.get_id()}, files={len(self.files)})'


    def to_ttl(self, f) -> None:
        dataset = self.get_dataset()
        project = dataset.get_project()

        for rawdata in self.files:
            f.write(f'bid:PRF{dataset.get_number()} jpost:hasRawData bid:{rawdata.get_id()} .\n')
            f.write(f'bid:{rawdata.get_id()} rdfs:label "{rawdata.get_file()}" ;\n')
            f.write(f'    foaf:page jrepo:{project.get_id()} ;\n')
            f.write(f'    a jpost:RawData .\n\n')

        

    @staticmethod
    def read_rawdata_list(dataset: DataSet, xml_path: str) -> RawDataList:
        tree = ET.parse(xml_path)
        root = tree.getroot()

        rawdata_list = RawDataList(dataset)

        files = root.findall('FileList/File')

        for file in files:
            name = None
            type = None
            
            name_tag = file.find('Name')
            if name_tag is not None:
                name = name_tag.text

            type_tag = file.find('Type')
            if type_tag is not None:
                type = type_tag.text

            if type is not None and type.lower() == 'raw':
                rawdata_list.add_file(name)

        return rawdata_list
