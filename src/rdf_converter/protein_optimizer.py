
from datetime import datetime
from pathlib import Path
import pathlib
import psutil
import datetime
import tqdm
from .utils.logging import get_logger
from .models.project import Project
from .models.dataset import DataSet
from .models.sample import Sample
from .models.fractionation import Fractionation
from .models.enzyme import Enzyme
from .models.msmode import MsMode
from .models.rawdata_list import RawDataList
from .models.peptide import Peptide
from .models.protein import Protein
from .models.isoform import Isoform
from .models.group import Group 
from .models.psm import Psm
from .models.spectrum import Spectrum

import requests

import os

logger = get_logger(__name__)


class ProteinOptimizer:
    def __init__(self, api_url: str, cache_dir: str):
        self.api_url = api_url
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(exist_ok=True)

    def list_datasets(self):
        url = f'{self.api_url}/dataset_id_list'
        response = requests.get(url)
        response.raise_for_status()
        return response.json()
    
    def clear_cache(self):
        for file in self.cache_dir.glob('*'):
            file.unlink()

    def save_cache(self, dataset_id: str, overwrite: bool = False) -> None:
        file_path = self.cache_dir / f'{dataset_id}.txt'

        if not file_path.exists() or overwrite:
            with open(file_path, 'w') as f:
                url = f'{self.api_url}/dataset_protein_pepseq_score_list?dataset_id={dataset_id}'
                response = requests.get(url)
                response.raise_for_status()
                objects = response.json()
                for object in objects:
                    dataset_id = object['dataset_id']
                    uniprot = object['uniprot']
                    pep_seq = object['pep_seq']
                    score = object['max_score']
                    f.write(f'{dataset_id}\t{uniprot}\t{pep_seq}\t{score}\n')


    def load_cache(self, dataset_ids: list[str], min_score: float) -> list[Protein]:
        protein_map = {}
        peptide_map = {}
        sequence_map = {}

        project = Project('9999')
        dataset = DataSet(project, '9999')

        proteins = []

        for dataset_id in dataset_ids:
            file_path = self.cache_dir / f'{dataset_id}.txt'
            if not file_path.exists():
                self.save_cache(dataset_id)

            with open(file_path, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) == 4:
                        dataset_id, uniprot, pep_seq, score = parts
                        score = int(score)

                        if score >= min_score:
                            protein = None
                            if uniprot in protein_map:
                                protein = protein_map[uniprot]
                            else:
                                protein = Protein(dataset, uniprot)
                                protein_map[uniprot] = protein
                                proteins.append(protein)

                            peptide = None
                            if pep_seq in peptide_map:
                                peptide = peptide_map[pep_seq]
                                if score > peptide.get_score():
                                    peptide.set_score(score)
                            else:
                                peptide = Peptide(dataset, pep_seq)
                                peptide_map[pep_seq] = peptide
                                peptide.set_score(score)

                            sequence_set = None
                            if uniprot in sequence_map:
                                sequence_set = sequence_map[uniprot]
                            else:
                                sequence_set = set()
                                sequence_map[uniprot] = sequence_set
                            if pep_seq not in sequence_set:
                                sequence_set.add(pep_seq)
                                protein.add_match(peptide, '', '')
        return proteins
    
    
    def optimize_proteins(self, dataset_ids: list[str], min_score: int, update_cache: bool) -> list[Protein]:
        logger.info('Creating cache')
        for i in tqdm.tqdm(range(len(dataset_ids))):
            dataset_id = dataset_ids[i]
            self.save_cache(dataset_id, update_cache)
        proteins = self.load_cache(dataset_ids, min_score)
        optimized_proteins = Protein.optimize(proteins)
        return optimized_proteins
    