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
class Peptide:
    dataset: DataSet
    id: str | None = None
    psms: list[Psm] = field(default_factory=list)
    sequence: str | None = None
    dummy: str | None = None
    unique: bool = False
    unique_at_mslevel: bool = False
    score: float = 0.0
    fdr: float = 0.0
    mod: str | None = None
    distinguishable_peptides: list[Peptide] = field(default_factory=list)

    counter: ClassVar[int] = 0


    def __init__(self, dataset: DataSet, sequence: str):
        Peptide.counter += 1
        self.dataset = dataset
        self.id = f'PEP{dataset.get_number()}_{Peptide.counter}'

        self.psms = []
        self.sequence = sequence
        self.set_dummy(sequence)
        self.distinguishable_peptides = []

        self.score = 0.0
        self.fdr = 0.0

    
    def get_dataset(self) -> DataSet:
        return self.dataset
    
    def get_id(self) -> str | None:
        return self.id
    
    def set_id(self, id: str) -> None:
        self.id = id
    
    def get_psms(self) -> list[Psm]:
        return self.psms
    
    def set_psms(self, psms: list[Psm]) -> None:
        self.psms = psms
    
    def get_sequence(self) -> str | None:
        return self.sequence
    
    def set_dummy(self, sequence: str) -> None:
        dummy = sequence.replace('I', 'J')
        dummy = dummy.replace('L', 'J')
        self.dummy = dummy

    def get_dummy(self) -> str | None:
        return self.dummy
    
    def is_unique(self) -> bool:
        return self.unique
    
    def set_unique(self, unique: bool) -> None:
        self.unique = unique

    def is_unique_at_mslevel(self) -> bool:
        return self.unique_at_mslevel   
    
    def set_unique_at_mslevel(self, unique_at_mslevel: bool) -> None:
        self.unique_at_mslevel = unique_at_mslevel

    def get_score(self) -> float:
        return self.score
    
    def set_score(self, score: float) -> None:
        self.score = score

    def get_fdr(self) -> float:
        return self.fdr
    
    def set_fdr(self, fdr: float) -> None:
        self.fdr = fdr

    def get_mod(self) -> str | None:
        return self.mod

    def set_mod(self, mod: str | None) -> None:
        self.mod = mod

    def get_distinguishable_peptides(self) -> list[Peptide]:
        return self.distinguishable_peptides

    def __str__(self):
        return f'Peptide(id={self.id}, sequence={self.sequence}, psms={len(self.psms)})'
    

    def to_ttl(self, f) -> None:
        f.write(f':{self.id} a jpost:Peptide ;\n')
        if self.is_unique_at_mslevel():
            f.write(f'    a jpost:UniquePeptideAtMsLevel ;\n')
        elif self.is_unique():
            f.write(f'    a jpost:SharedPeptideAtMsLevel ;\n')

        if self.is_unique():
            f.write(f'    a jpost:UniquePeptide ;\n')
        else:
            f.write(f'    a jpost:SharedPeptide ;\n')

        f.write('    sio:SIO_000216 [\n')
        f.write('        a jpost:UniScore ;\n')
        f.write(f'        sio:SIO_000300 {self.get_score()} ;\n')
        f.write('    ];\n')

        f.write('    sio:SIO_000216 [\n')
        f.write('        a obo:MS_1001364 ;\n')
        f.write(f'        sio:SIO_000300 {self.get_fdr()} ;\n')
        f.write('    ];\n')

        f.write('    jpost:hasSequence [\n')
        f.write('        a obo:MS_1001344 ;\n')
        f.write(f'        rdf:value "{self.get_sequence()}" ;\n')
        f.write('    ];\n')

        for psm in self.get_psms():
            f.write(f'    jpost:hasPsm :{psm.get_id()} ;\n')


        for indistinguishable in self.get_distinguishable_peptides():
            f.write(f'    jpost:hasIndistinguishablePeptide :{indistinguishable.get_id()} ;\n')
        f.write(f'    dct:identifier "{self.get_id()}" .\n\n')



    
    @staticmethod
    def read_peptides(dataset: DataSet, tsv_path: str) -> list[Peptide]:
        peptides: list[Peptide] = []

        if not os.path.exists(tsv_path):
            return peptides
        


        with open(tsv_path, 'r', encoding='utf-8') as f:
            reader = csv.reader(f, delimiter='\t')

            header = next(reader)
            header = [col.lower().replace(' ', '').replace('-', '').replace('_', '').strip() for col in header]

            index_hits = -1
            index_rts = -1
            index_scores = -1
            index_jpost_scores = -1
            index_masses = -1
            index_charges = -1
            index_files = -1
            index_scans = -1
            index_phospho_confirmed = -1
            index_sequence = -1
            index_mod = -1
            index_mod_detail = -1
            index_ambigious = -1
            index_fdr = -1
            index_title = -1
            index_calcmz = -1

            for i, col in enumerate(header):
                if col == 'hitpsmcount':
                    index_hits = i
                elif col == 'sameseqrtime':
                    index_rts = i
                elif col == 'sameseqplenghitscore':
                    index_scores = i
                elif col == 'sameseqjpostscore':
                    index_jpost_scores = i
                elif col == 'sameseqobsmass':
                    index_masses = i
                elif col == 'sameseqcharge':
                    index_charges = i
                elif col == 'sameseqrawfile':
                    index_files = i
                elif col == 'sameseqscanno':
                    index_scans = i
                elif col == 'sameseqphosphoconfimedsite':
                    index_phospho_confirmed = i
                elif col == 'seq':
                    index_sequence = i
                elif col == 'mod':
                    index_mod = i
                elif col == 'moddetail':
                    index_mod_detail = i
                elif col == 'sameseqphosphoambiguoussite':
                    index_ambigious = i
                elif col == 'pepfdr':
                    index_fdr = i
                elif col == 'title':
                    index_title = i
                elif col == 'calcmz':
                    index_calcmz = i

            for row in reader:
                hits = int(row[index_hits])
                rts = row[index_rts].split(',')
                scores = row[index_scores].split(',')
                jpost_scores = row[index_jpost_scores].split(',')
                masses = row[index_masses].split(',')
                charges = row[index_charges].split(',')
                files = row[index_files].split(',')
                scans = row[index_scans].split(',')
                phospho_confirmed = row[index_phospho_confirmed]
                sequence = row[index_sequence]
                mod = row[index_mod]
                mod_detail = row[index_mod_detail]
                title = row[index_title]
                calc_mz = float(row[index_calcmz])

                confirmed = []
                if phospho_confirmed is not None and phospho_confirmed != '':
                    confirmed = phospho_confirmed.split(',')

                max_jpost_score = 0.0

                phospho_ambigious = row[index_ambigious]
                ambigious = []
                if phospho_ambigious is not None and phospho_ambigious != '':
                    ambigious = phospho_ambigious.split(',')

                fdr = float(row[index_fdr])

                peptide = Peptide(dataset, sequence)
                peptide.set_mod(mod if mod.startswith('Label:') else '')

                psms = []
                for i in range(hits):
                    current_psms = []
                    if mod.startswith('Label:'):
                        psm1 = Psm(dataset)
                        psm1.set_mod('')
                        psm1.set_mod_detail('')
                        current_psms.append(psm1)

                        psm2 = Psm(dataset)
                        psm2.set_mod(mod)
                        psm2.set_mod_detail(mod_detail)
                        current_psms.append(psm2)
                    else:
                        psm = Psm(dataset)
                        psm.set_mod(mod)
                        psm.set_mod_detail(mod_detail)
                        current_psms.append(psm)

                    for psm in current_psms:
                        psm.set_sequence(sequence)
                        psm.set_title(title)
                        psm.set_calc_mz(calc_mz)
                        psm.set_fdr(fdr)

                        if i < len(rts):
                            psm.set_rt(rts[i])
                        if i < len(jpost_scores):
                            psm.set_jpost_score(jpost_scores[i])
                            try:
                                jpost_score = float(jpost_scores[i])
                                if jpost_score > max_jpost_score:
                                    max_jpost_score = jpost_score
                            except:
                                pass
                        if i < len(masses):
                            psm.set_obs_mz(masses[i])
                        if i < len(charges):
                            psm.set_charge(charges[i])
                        if i < len(files):
                            psm.set_raw_file(files[i])
                        if i < len(scans):
                            psm.set_scan(scans[i])
                        if i < len(scores):
                            tokens = scores[i].split('/')
                            if len(tokens) >= 4:
                                psm.get_score_map()['ev'] = tokens[2]
                                psm.get_score_map()[tokens[1]] = tokens[3]
                        if confirmed is not None and i < len(confirmed):
                            psm.set_phospho_confirmed(confirmed[i])
                        if ambigious is not None and i < len(ambigious):
                            tokens = ambigious[i].split('/')
                            for token in tokens:
                                if token != '':
                                    psm.get_phospho_ambiguous().append(token)
                        psm.set_peptide(peptide)
                        psms.append(psm)

                    
                peptide.set_score(max_jpost_score)
                peptide.set_fdr(fdr)
                peptide.set_psms(psms)
                peptides.append(peptide)

        return peptides
    
    @staticmethod
    def check_peptides(proteins: list[Protein], peptides: list[Peptide]) -> None:
        sequence_map = {}
        dummy_map = {}
        protein_map = {}
        count_map = {}

        for protein in proteins:
            protein_map[protein.get_uniprot()] = protein
            peptide_list = protein.search_peptides()
            sequence_set = set()
            dummy_set = set()
            for peptide in peptide_list:
                sequence_set.add(peptide.get_sequence())
                dummy_set.add(peptide.get_dummy())
                count_key = f'{peptide.get_dummy()}_{peptide.get_mod()}'
                if count_key not in count_map:
                    count_map[count_key] = [peptide]
                else:
                    count_map[count_key].append(peptide)

            sequence_map[protein.get_uniprot()] = sequence_set
            dummy_map[protein.get_uniprot()] = dummy_set

        for peptide in peptides:
            seq_count = 0
            dummy_count = 0
            for protein in proteins:
                seq_set = sequence_map.get(protein.get_uniprot(), set())
                dummy_set = dummy_map.get(protein.get_uniprot(), set())
                if peptide.get_sequence() in seq_set:
                    seq_count += 1
                if peptide.get_dummy() in dummy_set:
                    dummy_count += 1

            count_key = f'{peptide.get_dummy()}_{peptide.get_mod()}'
            if count_key in count_map:
                count_list = count_map.get(count_key)
                for p in count_list:
                    if p != peptide:
                        peptide.get_distinguishable_peptides().append(p)

            peptide.set_unique(seq_count == 1)
            peptide.set_unique_at_mslevel(dummy_count == 1)
