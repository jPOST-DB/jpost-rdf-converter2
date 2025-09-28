from __future__ import annotations

from dataclasses import dataclass, field
from turtle import up
from typing import Optional, List, Dict, Tuple
import logging
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .dataset import DataSet
    from .peptide import Peptide
    from .group import Group
from .isoform import Isoform

logger = logging.getLogger(__name__)



from dotenv import load_dotenv
import os
import subprocess
import pulp
import heapq


from pathlib import Path


@dataclass
class PeptideMatch:
    peptide: Peptide | None = None
    start: str | None = None
    end: str | None = None

    def __init__(self):
        self.peptide = None
        self.start = None
        self.end = None

    def set_peptide(self, peptide: Peptide) -> None:
        self.peptide = peptide

    def get_peptide(self) -> Peptide | None:
        return self.peptide
    
    def set_start(self, start: str) -> None:
        self.start = start

    def get_start(self) -> str | None:
        return self.start   
    
    def set_end(self, end: str) -> None:
        self.end = end

    def get_end(self) -> str | None:
        return self.end
    
    def to_ttl(self, f) -> None:
        f.write(f'    jpost:hasPeptideEvidence [\n')
        f.write(f'        a jpost:PeptideEvidence ;\n')
        f.write(f'        jpost:hasPeptide :{self.get_peptide().get_id()} ;\n')
        f.write(f'        faldo:location [\n')
        f.write(f'            a faldo:Region ;\n')
        f.write(f'            faldo:begin [\n')
        f.write(f'                a faldo:ExactPosition ;\n')
        f.write(f'                faldo:reference :{self.get_peptide().get_id()} ;\n')
        f.write(f'                faldo:position {self.get_start()} ;\n')
        f.write(f'            ] ;\n')
        f.write(f'            faldo:end [\n')
        f.write(f'                a faldo:ExactPosition ;\n')
        f.write(f'                faldo:reference :{self.get_peptide().get_id()} ;\n')
        f.write(f'                faldo:position {self.get_end()} ;\n')
        f.write(f'            ] ;\n')
        f.write(f'        ] ;\n')
        f.write(f'    ] ;\n')


@dataclass
class Protein:
    dataset: DataSet | None = None
    id: str | None = None
    uniprot: str | None = None
    peptides: List[PeptideMatch] = field(default_factory=list)
    type: str | None = None
    group: Group | None = None
    in_optimization_list: bool = False
    leading: bool = False
    anchor: bool = False
    same: bool = False
    subset: bool = False
    leading_protein: list[Protein] = field(default_factory=list)
    isoforms: list[Isoform] = field(default_factory=list)

    def __init__(self, dataset: DataSet, uniprot: str):
        self.dataset = dataset
        self.id = f'PRT{dataset.get_number()}_{uniprot}'
        self.uniprot = uniprot
        self.peptide_matches = []
        self.type = None
        self.group = None
        self.in_optimization_list = False
        self.leading = False
        self.anchor = False
        self.same = False
        self.subset = False
        self.leading_protein = []
        self.isoforms = []

    def get_dataset(self) -> DataSet | None:
        return self.dataset
    
    def get_id(self) -> str | None:
        return self.id
    
    def get_uniprot(self) -> str | None:
        return self.uniprot
    
    def get_peptide_matches(self) -> list[PeptideMatch]:
        return self.peptide_matches
    
    def get_type(self) -> str | None:
        return self.type

    def set_type(self, type: str) -> None:
        self.type = type

    def get_group(self) -> Group | None:
        return self.group
    
    def set_group(self, group: Group) -> None:
        self.group = group

    def add_match(self, peptide: Peptide, start: str, end: str) -> None:
        match = PeptideMatch()
        match.set_peptide(peptide)
        match.set_start(start)
        match.set_end(end)
        self.peptide_matches.append(match)

    def is_in_optimization_list(self) -> bool:
        return self.in_optimization_list
    
    def set_in_optimization_list(self, in_optimization_list: bool) -> None:
        self.in_optimization_list = in_optimization_list

    def is_leading(self) -> bool:
        return self.leading
    
    def set_leading(self, leading: bool) -> None:
        self.leading = leading

    def is_anchor(self) -> bool:
        return self.anchor
    
    def set_anchor(self, anchor: bool) -> None:
        self.anchor = anchor

    def is_same(self) -> bool:
        return self.same
    
    def set_same(self, same: bool) -> None:
        self.same = same

    def is_subset(self) -> bool:
        return self.subset
    
    def set_subset(self, subset: bool) -> None:
        self.subset = subset

    def get_leading_proteins(self) -> list[Protein]:
        return self.leading_protein
    
    def set_isoform(self, isoform: list[Isoform]) -> None:
        self.isoform = isoform

    def get_isoforms(self) -> list[Isoform]:
        return self.isoforms
    
    def is_included(self, protein: Protein) -> bool:
        sequence_set = set()
        for match in protein.get_peptide_matches():
            sequence_set.add(match.get_peptide().get_dummy())

        found = False
        for match in self.get_peptide_matches():
            if match.get_peptide().get_dummy() in sequence_set:
                found = True

        return found
    
    def search_peptides(self) -> list[Peptide]:
        peptides = []
        sequence_set = set()

        for match in self.get_peptide_matches():
            sequence = match.get_peptide().get_dummy()
            if sequence not in sequence_set:
                sequence_set.add(sequence)
                peptides.append(match.get_peptide())

        isoform = self.get_isoforms()
        for iso in isoform:
            for match in iso.get_peptide_matches():  
                sequence = match.get_peptide().get_dummy()
                if sequence not in sequence_set:
                    sequence_set.add(sequence)
                    peptides.append(match.get_peptide())

        return peptides


    def __str__(self):
        return f'Protein(id={self.id}, uniprot={self.uniprot}, peptides={len(self.peptide_matches)})'
    
    def to_ttl(self, f) -> None:
        f.write(f':{self.get_id()} a jpost:Protein ;\n')
        f.write(f'    dct:identifier "{self.get_id()}" ;\n')
        if self.get_group() is not None:
            f.write(f'    jpost:inProteinGroup bid:{self.get_group().get_id()} ;\n')
        f.write(f'    rdfs:label "{self.get_uniprot()}" ;\n')
        f.write(f'    rdfs:seeAlso idup:{self.get_uniprot()} ;\n')
        f.write(f'    rdfs:seeAlso uniprot:{self.get_uniprot()} ;\n')
        f.write(f'    jpost:hasDatabaseSequence uniprot:{self.get_uniprot()} ;\n')

        if self.is_in_optimization_list():
            f.write(f'    a jpost:RepresentativeIsoform ;\n')
        if self.is_anchor():
            f.write(f'    a obo:MS_1001591 ;\n')
        if self.is_leading():
            f.write(f'    a obo:MS_1002401 ;\n')
        elif self.is_same():
            f.write(f'    a obo:MS_1001595 ;\n')
        elif self.is_subset():
            f.write(f'    a obo:MS_1001597 ;\n')
        else:
            f.write(f'    a obo:MS_1001599 ;\n')

        for leading in self.get_leading_proteins():
            f.write(f'    jpost:hasLeadingProtein :{leading.get_id()} ;\n')

        f.write(f'    sio:SIO_000216 [\n')
        f.write(f'        a obo:MS_1001097 ;\n')
        f.write(f'        sio:SIO_000300 :{len(self.get_peptide_matches())} ;\n')
        f.write(f'    ] ;\n')

        f.write(f'    sio:SIO_000216 [\n')
        f.write(f'        a obo:MS_1002153 ;\n')
        f.write(f'        sio:SIO_000300 :{len(self.get_peptide_matches())} ;\n')
        f.write(f'    ] ;\n')        

        for match in self.get_peptide_matches():
            match.to_ttl(f)
        f.write(f'    a jpost:Protein .\n')
        f.write(f'\n')


    @staticmethod
    def create_db_index(fasta_path: str, work_dir: str) -> Path:
        load_dotenv()
        java_bin = os.getenv('JAVA_BIN', 'java')
        peptide_match_jar = os.getenv('PEPTIDEMATCH_JAR', './lib/PeptideMatchCMD.jar')  
        
        dir = Path(work_dir)
        db_index = dir / 'db_index'

        args = ['-a', 'index', '-d', fasta_path, '-i', db_index.resolve()]

        cmd = [java_bin, '-jar', peptide_match_jar] + args

        logger.info(f'Creating PeptideMatch DB index: {cmd}')

        with open(dir / 'db_index.log', 'w') as log_file:
            subprocess.run(cmd, stdout=log_file, stderr=subprocess.STDOUT, check=True)
    
        return db_index.resolve()
    

    @staticmethod
    def write_sequences_file(peptides: list[Peptide], work_dir: str) -> Path:
        dir = Path(work_dir)
        peptide_file = dir / 'peptides.txt'
        with open(peptide_file, 'w') as f:
            for peptide in peptides:
                f.write(f'{peptide.get_sequence()}\n')
        return peptide_file.resolve()
    

    @staticmethod
    def execute_peptide_match(peptides: list[Peptide], db_index: Path, work_dir: str) -> Path:
        load_dotenv()

        java_bin = os.getenv('JAVA_BIN', 'java')
        peptide_match_jar = os.getenv('PEPTIDEMATCH_JAR', './lib/PeptideMatchCMD.jar')  
        
        dir = Path(work_dir)
        output_path = dir / 'peptide_matches.txt'

        peptide_list_file = Protein.write_sequences_file(peptides, work_dir)

        args = [
            '-a',
            'query',
            '-l',
            '-Q',
            str(peptide_list_file),
            '-i',
            str(db_index),
            '-o',
            str(output_path)
        ]

        cmd = [java_bin, '-jar', peptide_match_jar] + args

        logger.info(f'Executing PeptideMatch: {cmd}')

        with open(dir / 'peptide_match.log', 'w') as log_file:
            subprocess.run(cmd, stdout=log_file, stderr=subprocess.STDOUT, check=True)
    
        return output_path.resolve()



    
    @staticmethod
    def get_protein_list(peptides: list[Peptide], work_dir: str, fasta_path: str) -> list[Protein]:
        proteins = []

        db_index = Protein.create_db_index(fasta_path, work_dir)
        peptide_match_file = Protein.execute_peptide_match(peptides, db_index, work_dir)

        peptide_map = {peptide.get_sequence(): peptide for peptide in peptides}
        protein_map = {}
        sequence_set = set(peptide_map.keys())

        with open(peptide_match_file, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    parts = line.strip().split('\t')
                    if len(parts) >= 5:
                        sequence = parts[0]
                        uniprot = parts[1]
                        uniprot_id = uniprot.split('|')[1] if '|' in uniprot else uniprot
                        start = parts[2]
                        end = parts[3]

                        peptide = peptide_map.get(sequence)
                        protein = protein_map.get(uniprot_id)
                        if protein is None:
                            protein = Protein(peptide.get_dataset(), uniprot_id)
                            protein_map[uniprot_id] = protein
                            proteins.append(protein)

                        protein.add_match(peptide, start, end)
                        if sequence in sequence_set:
                            sequence_set.remove(sequence)

        last_peptides = [pep for pep in peptides if not pep.get_sequence() in sequence_set]

        counter = 0
        for peptide in last_peptides:
            counter += 1
            id = f'PEP{peptide.get_dataset().get_number()}_{counter}'
            peptide.set_id(id)

        return {
            'peptides': last_peptides,
            'proteins': proteins
        }
    
    @staticmethod
    def create_isoforms(proteins: list[Protein], optimized_proteins: list[Protein]) -> list[Protein]:
        new_proteins = []
        isoforms = []

        optimized_set = set()
        for protein in optimized_proteins:
            optimized_set.add(protein.get_uniprot())

        protein_map = {protein.get_uniprot(): protein for protein in proteins if protein.get_uniprot().find('-') < 0}

        for protein in optimized_proteins:
            uniprot = protein.get_uniprot()
            index = uniprot.rfind('-')
            if index < 0:
                base_protein = protein
                base_protein.set_in_optimization_list(uniprot in optimized_set)
                new_proteins.append(base_protein)
                protein_map[base_protein.get_uniprot()] = base_protein
            else:
                uniprot = uniprot[:index]
                base_protein = protein_map.get(uniprot)
                if base_protein is not None:
                    base_protein = Protein(protein.get_dataset(), uniprot)
                    protein_map[uniprot] = base_protein
                    new_proteins.append(base_protein)

                    isoform = Isoform()
                    isoform.set_protein(base_protein)
                    base_protein.get_isoforms().append(isoform)
                    isoform.set_uniprot(protein.get_uniprot())
                    isoform.set_id(protein.get_id().replace('PRT', 'ISO'))
                    isoform.set_in_optimization_list(protein.get_uniprot() in optimized_set)
                    for match in protein.get_peptide_matches():
                        isoform.get_peptide_matches().append(match)
                    isoforms.append(isoform)

        return new_proteins, isoforms
    

    @staticmethod
    def check_proteins(proteins: list[Protein]) -> None:
        protein_map = {}
        sequence_map = {}

        for protein in proteins:
            uniprot = protein.get_uniprot()
            protein_map[uniprot] = protein

            leading = protein.is_in_optimization_list()
            for isoform in protein.get_isoforms():
                if isoform.is_in_optimization_list():
                    leading = True

            protein.set_leading(leading)
            protein.set_anchor(False)
            protein.set_same(False)
            protein.set_subset(False)

            if leading:
                peptides = protein.search_peptides()
                sequence_map[uniprot] = peptides

        for protein in proteins:
            if not protein.is_leading():
                peptides = protein.search_peptides()

                for uniprot in protein_map.keys():
                    leading_peptides = sequence_map[uniprot]
                    subset = True
                    count = 0
                    for peptide in peptides:
                        for leading_peptide in leading_peptides:
                            if peptide.get_dummy() == leading_peptide.get_dummy():
                                count += 1
                            else:
                                subset = False
                    
                    if subset:
                        if count == len(leading_peptides):
                            protein.set_same(True)
                            protein_map[uniprot].set_anchor(True)
                        else:
                            protein.set_subset(True)
                            protein.get_leading_proteins().append(protein_map[uniprot])


    @staticmethod    
    def optimize(proteins: list[Protein]) -> list[Protein]:
        load_dotenv()
        protein_parameter = int(os.getenv('PROTEIN_PARAMETER', '5000'))
        peptide_parameter = int(os.getenv('PEPTIDE_PARAMETER', '10000'))    

        peptide_count = 0
        for protein in proteins:
            peptide_count += len(protein.get_peptide_matches())
    
        optimized_proteins = None
        if len(proteins) > protein_parameter or peptide_count > peptide_parameter:
            optimized_proteins = Protein.solve_set_cover_by_greedy(proteins)
        else:
            optimized_proteins = Protein.solve_set_cover_by_ilp(proteins)

        return optimized_proteins


    @staticmethod
    def solve_set_cover_by_greedy(protein_list: list[Protein]) -> list[Protein]:
        protein_to_peptides: dict[str, set] = {}
        protein_lookup: dict[str, "Protein"] = {}
        all_peptides: set = set()

        for protein in protein_list:
            uniprot_id = protein.get_uniprot()
            peptide_set = {match.get_peptide().get_dummy() for match in protein.get_peptide_matches()}
            if len(peptide_set) > 0:
                protein_to_peptides[uniprot_id] = peptide_set
                protein_lookup[uniprot_id] = protein
                all_peptides |= peptide_set

        uncovered_peptides = set(all_peptides)
        if len(uncovered_peptides) == 0:
            return [], 0.0

        heap: list[tuple[int, int, str]] = []  # (-カバー数, 世代番号, uniprot_id)
        current_generation = 0
        coverage_count_cache: dict[str, int] = {}

        for uniprot_id, peptide_set in protein_to_peptides.items():
            coverage_count = len(peptide_set)
            if coverage_count > 0:
                coverage_count_cache[uniprot_id] = coverage_count
                heapq.heappush(heap, (-coverage_count, current_generation, uniprot_id))

        selected_uniprots: list[str] = []
        total_cost = 0.0

        while len(uncovered_peptides) > 0 and len(heap) > 0:
            neg_coverage, generation, candidate_uniprot = heapq.heappop(heap)

            recalculated_coverage = len(protein_to_peptides[candidate_uniprot] & uncovered_peptides)
            is_generation_outdated = generation != current_generation
            is_candidate_useful = recalculated_coverage > 0

            if is_generation_outdated and is_candidate_useful:
                coverage_count_cache[candidate_uniprot] = recalculated_coverage
                heapq.heappush(heap, (-recalculated_coverage, current_generation, candidate_uniprot))
            elif is_candidate_useful:
                selected_uniprots.append(candidate_uniprot)
                total_cost += 1.0
                newly_covered_peptides = protein_to_peptides[candidate_uniprot] & uncovered_peptides
                uncovered_peptides -= newly_covered_peptides
                current_generation += 1

        selected_proteins = [protein_lookup[u] for u in selected_uniprots]
        return selected_proteins

    @staticmethod
    def solve_set_cover_by_ilp(protein_list: list[Protein]) -> list[Protein]:
        protein_to_peptides: dict[str, set] = {}
        protein_lookup: dict[str, "Protein"] = {}
        all_peptides: set = set()

        for protein in protein_list:
            uniprot_id = protein.get_uniprot()
            peptide_set = {match.get_peptide().get_dummy() for match in protein.get_peptide_matches()}
            if len(peptide_set) > 0:
                protein_to_peptides[uniprot_id] = peptide_set
                protein_lookup[uniprot_id] = protein
                all_peptides |= peptide_set

        if len(all_peptides) == 0:
            return []

        prob = pulp.LpProblem("Set_Cover_Problem", pulp.LpMinimize)

        protein_vars: dict[str, pulp.LpVariable] = {}
        for uniprot_id in protein_to_peptides.keys():
            var = pulp.LpVariable(f'P_{uniprot_id}', cat='Binary')
            protein_vars[uniprot_id] = var

        prob += pulp.lpSum(protein_vars.values()), "Total_Proteins_Selected"

        for peptide in all_peptides:
            prob += (
                pulp.lpSum(
                    protein_vars[uniprot_id]
                    for uniprot_id, peptides in protein_to_peptides.items()
                    if peptide in peptides
                ) >= 1,
                f"Cover_Peptide_{peptide}"
            )

        prob.solve(pulp.PULP_CBC_CMD(msg=False))

        selected_proteins = [
            protein_lookup[uniprot_id]
            for uniprot_id, var in protein_vars.items()
            if var.varValue == 1.0
        ]

        return selected_proteins