from datetime import datetime
from pathlib import Path
import pathlib
import psutil
import datetime
import logging
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
from .models.pep import Pep

import os

logger = get_logger(__name__)

class ProjectConverter:
    def __init__(self, rev: str, meta_data_path: str, out_path: str):
        self.rev = rev
        self.meta_data_path = meta_data_path
        self.out_path = out_path

    def convert(self):
        logger.info(f'Starting project conversion for rev: {self.rev}')
        project = Project.read_project(self.meta_data_path)
        project.set_id(self.rev)

        self.write_ttl(self.out_path, project)
        

    def write_ttl(self, out_path: str, project: Project):
        with open(out_path, 'w', encoding='utf-8') as f:
            self.write_header(f)
            project.to_ttl(f)


    def write_header(self, f) -> None:
        headers = [
            '@prefix rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> .',
			'@prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#> .', '@prefix dct: <http://purl.org/dc/terms/> .',
			'@prefix uniprot: <http://purl.uniprot.org/uniprot/> .',
			'@prefix isoforms: <http://purl.uniprot.org/isoforms/> .',
			'@prefix idup: <http://identifiers.org/uniprot/> .',
			'@prefix taxonomy: <http://identifiers.org/taxonomy/> .',
			'@prefix obo: <http://purl.obolibrary.org/obo/> .',
			'@prefix ncit: <http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#> .',
			'@prefix unimod: <http://www.unimod.org/obo/unimod.obo#> .',
			'@prefix sio: <http://semanticscience.org/resource/> .', '@prefix foaf: <http://xmlns.com/foaf/0.1/> .',
			'@prefix faldo: <http://biohackathon.org/resource/faldo#> .',
			'@prefix skos: <http://www.w3.org/2004/02/skos/core#> .',
			'@prefix px: <https://github.com/PX-RDF/ontology/blob/master/px.owl#> .',
			'@prefix pxd: <http://proteomecentral.proteomexchange.org/dataset/> .',
			'@prefix jpost: <http://rdf.jpostdb.org/ontology/jpost.owl#> .',
			'@prefix jrepo: <https://repository.jpostdb.org/entry/> .',
			'@prefix bid: <http://rdf.jpostdb.org/bid/> .',
			'@prefix vcard: <http://www.w3.org/2006/vcard/ns#> .',
			'@prefix xsd: <http://www.w3.org/2001/XMLSchema#> .',
			'@prefix : <http://rdf.jpostdb.org/entry/> .'
        ]

        for line in headers:
            f.write(f'{line}\n')
        f.write('\n')

