from __future__ import annotations

from dataclasses import dataclass, field
from ..utils.string_tool import is_not_empty
import logging

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .dataset import DataSet


logger = logging.getLogger(__name__)

@dataclass
class Fasta:
    title: str
    sequence: str

    def __init(self, title: str, sequence: str):
        self.title = title
        self.sequence = sequence

    def get_title(self) -> str:
        return self.title
    
    def get_sequence(self) -> str:
        return self.sequence
    
    @staticmethod
    def read_fasta(fasta_path: str) -> list[Fasta]:
        fastas = []
        with open(fasta_path, 'r') as fr:
            title = ''
            seq_lines = []
            for line in fr:
                line = line.strip()
                if line.startswith('>'):
                    if is_not_empty(title) and len(seq_lines) > 0:
                        sequence = ''.join(seq_lines)
                        fasta = Fasta(title, sequence)
                        fastas.append(fasta)
                    title = line[1:].strip()
                    seq_lines = []
                else:
                    seq_lines.append(line.strip())
            if is_not_empty(title) and len(seq_lines) > 0:
                sequence = ''.join(seq_lines)
                fasta = Fasta(title, sequence)
                fastas.append(fasta)
        return fastas
    
