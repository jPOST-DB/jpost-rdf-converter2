from __future__ import annotations

from dataclasses import dataclass, field
import logging

from statistics import mode
from xml.etree import ElementTree as ET

from typing import TYPE_CHECKING

from ..utils.string_tool import is_not_empty

if TYPE_CHECKING:
    from .dataset import DataSet

logger = logging.getLogger(__name__)

@dataclass
class MsMode:
    dataset: DataSet
    id: str | None = None
    instrument: str | None = None
    instrument_mode: str | None = None
    quantification_method: str | None = None
    purpose: str | None = None
    comments: str | None = None

    def __init__(self, dataset: DataSet):
        self.dataset = dataset
        self.dataset.set_ms_mode(self)
        self.id = f'MSM{dataset.get_number()}'

    def get_dataset(self) -> DataSet:
        return self.dataset
    
    def get_id(self) -> str | None: 
        return self.id
    
    def get_instrument(self) -> str | None:
        return self.instrument
    
    def set_instrument(self, instrument: str) -> None:
        self.instrument = instrument

    def get_instrument_mode(self) -> str | None:
        return self.instrument_mode
    
    def set_instrument_mode(self, instrument_mode: str) -> None:
        self.instrument_mode = instrument_mode

    def get_quantification_method(self) -> str | None:
        return self.quantification_method
    
    def set_quantification_method(self, quantification_method: str) -> None:
        self.quantification_method = quantification_method

    def get_purpose(self) -> str | None:
        return self.purpose
    
    def set_purpose(self, purpose: str) -> None:
        self.purpose = purpose

    def get_comments(self) -> str | None:
        return self.comments
    
    def set_comments(self, comments: str) -> None:
        self.comments = comments

    def __str__(self):
        return f'MsMode(id={self.id}, instrument={self.instrument}, instrument_mode={self.instrument_mode}, quantification_method={self.quantification_method}, purpose={self.purpose}, comments={self.comments})'
    
    def to_ttl(self, f) -> None:
        ms_mode = self
        dataset = self.get_dataset()
        
        f.write(f'bid:PRF{dataset.get_number()} jpost:hasMsMode bid:{ms_mode.get_id()} .\n')
        f.write(f'bid:{ms_mode.get_id()}\n')

        instrument = ms_mode.get_instrument()
        if is_not_empty(instrument):
            f.write(f'    jpost:instrument obo:{instrument.replace(":", "_")} ;\n')

        instrument_mode = ms_mode.get_instrument_mode()
        if is_not_empty(instrument_mode):
            f.write(f'    jpost:instrumentMode jpost:{instrument_mode.replace(":", "_")} ;\n')

        purpose = ms_mode.get_purpose()
        if is_not_empty(purpose):
            f.write(f'    jpost:purpose jpost:{purpose.replace(":", "_")} ;\n')

        quantification_mode = ms_mode.get_quantification_method()
        if is_not_empty(quantification_mode):
            f.write(f'    jpost:quantificationMethod obo:{quantification_mode.replace(":", "_")} ;\n')

        comment = ms_mode.get_comments()
        if is_not_empty(comment):
            f.write(f'    rdfs:comment "{comment.replace(":", "_")}" ;\n')  

        f.write(f'    a jpost:MsMode .\n\n')


    @staticmethod
    def read_ms_mode(dataset: DataSet, xml_path: str) -> MsMode:
        ms_mode = MsMode(dataset)

        tree = ET.parse(xml_path)
        root = tree.getroot()

        tag = root.find('FileList/File/Profile/MS_mode')
        if tag is not None:
            instrument = tag.find('instrument')
            if instrument is not None:
                ms_mode.set_instrument(instrument.get('id'))
            
            instrument_mode = tag.find('instrumentMode')
            if instrument_mode is not None:
                ms_mode.set_instrument_mode(instrument_mode.get('id'))

            purpose = tag.find('purpose')
            if purpose is not None:
                ms_mode.set_purpose(purpose.get('id'))

            note = tag.find('note')
            if note is not None and note.text is not None:
                text = note.text
                lines = text.split('\n')
                for line in lines:
                    tokens = line.split('|')
                    if len(tokens) >= 2:
                        key = tokens[0].strip()
                        if key == 'Quantification Method':
                            ms_mode.set_quantification_method(tokens[2].strip())
                        elif key == 'Note':
                            ms_mode.set_comments(tokens[1].strip())

        return ms_mode


