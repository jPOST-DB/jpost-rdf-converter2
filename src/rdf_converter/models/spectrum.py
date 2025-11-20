from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional, List

import logging

from ..utils.string_tool import is_not_empty

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .dataset import DataSet
    from .rawdata_list import RawData

logger = logging.getLogger(__name__)

@dataclass
class Spectrum:
    rawdata: RawData
    scan: str
    id: str | None = None

    def __init__(self, rawdata: RawData, scan: str):
        self.rawdata = rawdata
        self.scan = scan
        self.id = f'{rawdata.get_id().replace("RAW", "SPC")}_{scan}'

    def get_rawdata(self) -> RawData:
        return self.rawdata
    
    def get_scan(self) -> str:
        return self.scan
    
    def get_id(self) -> str | None:
        return self.id

    def to_ttl(self, f) -> None:
        f.write(f'bid:{self.id} a jpost:Spectrum ;\n')
        f.write(f'    rdfs:label "{self.id}" ;\n')
        f.write(f'    jpost:inRawData bid:{self.rawdata.get_id()} .\n\n')