from __future__ import annotations

from dataclasses import dataclass, field
from typing import ClassVar

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from .dataset import DataSet
    from .protein import Protein

from .psm import Psm

import os
import csv

@dataclass
class Pep :
    dataset: DataSet
    uniscore: int | None = None    
    hit_count: int | None = None
    decoy_hit_count: int | None = None

    def __init__(self, dataset, uniscore: int):
        self.dataset = dataset
        self.uniscore = uniscore

    def get_dataset(self) -> DataSet:
        return self.dataset
    
    def get_uniscore(self) -> int | None:
        return self.uniscore
    
    def set_hit_count(self, count: int) -> None:
        self.hit_count = count

    def get_hit_count(self) -> int | None:
        return self.hit_count
    
    def set_decoy_hit_count(self, count: int) -> None:
        self.decoy_hit_count = count

    def get_decoy_hit_count(self) -> int | None:
        return self.decoy_hit_count
    
    def to_ttl(self, f) -> None:
        f.write(f':{self.get_dataset().get_id()} jpost:hasPEP [\n')
        f.write(f'    sio:SIO_000216 [\n')
        f.write(f'        a jpost:UniScore ;\n')
        f.write(f'        sio:SIO_000300 {self.get_uniscore()} ;\n')
        f.write(f'    ] ;\n')
        if self.get_hit_count() is not None:
            f.write(f'    sio:SIO_000216 [\n')
            f.write(f'        a jpost:NormalHitCount ;\n')
            f.write(f'        sio:SIO_000300 {self.get_hit_count()} ;\n')
            f.write(f'    ] ;\n')
        if self.get_decoy_hit_count() is not None:
            f.write(f'    sio:SIO_000216 [\n')
            f.write(f'        a jpost:DecoyHitCount ;\n')
            f.write(f'        sio:SIO_000300 {self.get_decoy_hit_count()} ;\n')
            f.write(f'    ] ;\n')
        f.write(f'] .\n')
    
    @staticmethod
    def read_pep(dataset: DataSet, pep_path: str) -> list[Pep] | None:
        peps = []
        with open(pep_path, 'r') as fr:
            reader = csv.reader(fr, delimiter='\t')
            headers = next(reader)
            uniscore_idx = headers.index('jPostScore') if 'jPostScore' in headers else -1
            hit_count_idx = headers.index('NormalHitCount') if 'NormalHitCount' in headers else -1
            decoy_hit_count_idx = headers.index('DecoyHitCount') if 'DecoyHitCount' in headers else -1

            print(f'UniScore index: {uniscore_idx}, HitCount index: {hit_count_idx}, DecoyHitCount index: {decoy_hit_count_idx}')

            pep = None
            for row in reader:
                if uniscore_idx >= 0:
                    uniscore = float(row[uniscore_idx])
                    if uniscore.is_integer():
                        uniscore = int(uniscore)
                        pep = Pep(dataset, uniscore)
                        if hit_count_idx >= 0:
                            pep.set_hit_count(int(float(row[hit_count_idx])))
                        if decoy_hit_count_idx >= 0:
                            pep.set_decoy_hit_count(int(float(row[decoy_hit_count_idx])))
                        peps.append(pep)
        return peps
    
    