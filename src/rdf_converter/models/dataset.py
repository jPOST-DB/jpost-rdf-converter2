from __future__ import annotations

from dataclasses import dataclass, field
import logging

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .project import Project
    from .sample import Sample
    from .fractionation import Fractionation
    from .enzyme import Enzyme
    from .msmode import MsMode
    from .rawdata_list import RawDataList

from ..utils.string_tool import is_not_empty    

logger = logging.getLogger(__name__)

@dataclass
class DataSet:
    project: Project
    id: str | None = None
    number: str | None = None
    sample: Sample | None = None
    fractionation: Fractionation | None = None
    enzyme: Enzyme | None = None
    ms_mode: MsMode | None = None
    rawdata_list: RawDataList | None = None

    def __init__(self, project: Project, branch: str):
        self.project = project
        project.set_dataset(self)

        self.number = f'{project.get_project_number()}_{branch}'
        self.id = f'DS{self.number}'


    def get_project(self) -> Project:
        return self.project
    
    def get_id(self) -> str | None:
        return self.id
    
    def get_sample(self) -> Sample | None:
        return self.sample
    
    def set_sample(self, sample: Sample) -> None:
        self.sample = sample

    def get_fractionation(self) -> Fractionation | None:
        return self.fractionation
    
    def set_fractionation(self, fractionation: Fractionation) -> None:
        self.fractionation = fractionation

    def get_enzyme(self) -> Enzyme | None:
        return self.enzyme
    
    def set_enzyme(self, enzyme: Enzyme) -> None:
        self.enzyme = enzyme

    def get_ms_mode(self) -> MsMode | None:
        return self.ms_mode
    
    def set_ms_mode(self, ms_mode: MsMode) -> None:
        self.ms_mode = ms_mode

    def get_rawdata_list(self) -> RawDataList | None:
        return self.rawdata_list
    
    def set_rawdata_list(self, rawdata_list: RawDataList) -> None:
        self.rawdata_list = rawdata_list

    def get_number(self) -> str | None:
        return self.number
    

    def to_ttl(self, f) -> None:
        project = self.get_project()
        f.write(f':{project.get_id()} jpost:hasDataset :{self.get_id()} .\n\n')
        f.write(f':{self.get_id()} a jpost:DataSet ;\n')
        f.write(f'    jpost:hasProfile bid:PRF{self.get_number()} ;\n')
        f.write(f'    dct:identifier "{self.get_id()}" .\n\n')
        f.write(f'bid:PRF{self.get_number()} a jpost:Profile .\n\n')

        if self.sample is not None:
            self.sample.to_ttl(f)

        if self.fractionation is not None:
            self.fractionation.to_ttl(f)

        if self.enzyme is not None:
            self.enzyme.to_ttl(f)

        if self.ms_mode is not None:
            self.ms_mode.to_ttl(f)

        if self.rawdata_list is not None:
            self.rawdata_list.to_ttl(f)



