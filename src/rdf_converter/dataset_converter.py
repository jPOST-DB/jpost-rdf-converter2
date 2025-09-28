from datetime import datetime
from pathlib import Path
import pathlib
import psutil
import datetime
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

import os

logger = get_logger(__name__)

class DatasetConverter:
    def __init__(
            self, 
            rev: str,
            branch: str,
            tsv_path: str, 
            fasta_path: str,
            meta_path: str,
            result_dir: str,
            ttl_path: str,
            peptidematch_jar: str,
            java_bin: str
    ):
        self.rev = rev
        self.branch = branch
        self.tsv_path = Path(tsv_path)
        self.fasta_path = Path(fasta_path)
        self.meta_path = Path(meta_path)
        self.result_dir = Path(result_dir)
        self.ttl_path = Path(ttl_path)
        self.peptidematch_jar = Path(peptidematch_jar)
        self.java_bin = java_bin


    def get_work_folder(self) -> pathlib.Path:
        tmp_dir = pathlib.Path('tmp')
        tmp_dir.mkdir(exist_ok=True)

        pid = str(psutil.Process().pid)

        now = datetime.datetime.now()
        folder_name = now.strftime(f'%Y%m%d%H%M%S%f')[:-3] + f'_{pid}'

        folder_path = tmp_dir / folder_name
        folder_path.mkdir(exist_ok=True)

        return folder_path
    

    def convert(self) -> None:
        work_dir = self.get_work_folder().resolve()
        logger.info(f'Working directory: {work_dir}')

        project = Project.read_project(str(self.meta_path))
        project.set_id(self.rev)
        dataset = DataSet(project, self.branch)

        sample = Sample.read_sample(dataset, str(self.meta_path))
        logger.info(f'{sample}')

        fractionation = Fractionation.read_fractionation(dataset, str(self.meta_path))
        logger.info(f'{fractionation}')
        
        enzyme = Enzyme.read_enzyme(dataset, str(self.meta_path))
        logger.info(f'{enzyme}')

        ms_mode = MsMode.read_ms_mode(dataset, str(self.meta_path))
        logger.info(f'{ms_mode}')

        raw_data_list = RawDataList.read_rawdata_list(dataset, str(self.meta_path))
        logger.info(f'{raw_data_list}')

        peptides = Peptide.read_peptides(dataset, str(self.tsv_path))        
        logger.info(f'Peptides: {len(peptides)}')

        psms, spectra = Psm.get_psms(peptides)
        logger.info(f'PSMs: {len(psms)}')
        logger.info(f'Spectra: {len(spectra)}')

        protein_pair = Protein.get_protein_list(peptides, str(work_dir), str(self.fasta_path))
        proteins = protein_pair['proteins']
        peptides = protein_pair['peptides']
        logger.info(f'Hit Proteins: {len(proteins)}, Hit Peptides: {len(peptides)}')

        optimized_proteins = Protein.optimize(proteins)
        optimization_file = work_dir / 'optimization.txt'
        with open(optimization_file, 'w') as f:
            for protein in optimized_proteins:
                f.write(f'{protein.get_uniprot()}\n')
        logger.info(f'Optimized Proteins: {len(optimized_proteins)}')

        proteins, isoforms = Protein.create_isoforms(proteins, optimized_proteins)
        logger.info(f'Proteins: {len(proteins)}, Isoforms: {len(isoforms)}')

        groups = Group.create_groups(dataset, proteins, optimized_proteins)
        logger.info(f'Groups: {len(groups)}')

        Protein.check_proteins(proteins)
        Peptide.check_peptides(proteins, peptides)

        self.write_ttl(self.ttl_path, project, dataset, peptides, proteins, optimized_proteins, isoforms, groups, psms, spectra)


    def write_ttl(
            self, 
            ttl_path: pathlib.Path, 
            project: Project,
            dataset: DataSet,
            peptides: list[Peptide],
            proteins: list[Protein],
            optimized_proteins: list[Protein],
            isoforms: list[Isoform],
            groups: list[Group],
            psms: list[Psm],
            spectra: list[Spectrum]
    ) -> None:
        with open(ttl_path, 'w', encoding='utf-8') as f:
            self.write_header(f)
            project.to_ttl(f)
            dataset.to_ttl(f)

            for spectrum in spectra:
                spectrum.to_ttl(f)

            all_not_found = []
            for psm in psms:
                not_found = psm.to_ttl(f)
                all_not_found.extend(not_found)

            all_not_found.sort()
            all_not_found = list(dict.fromkeys(all_not_found))
            for not_found_modification in all_not_found:
                logger.warning(f'Modification not found: {not_found_modification}') 

            for peptide in peptides:
                peptide.to_ttl(f)

            for group in groups:
                group.to_ttl(f)

            for protein in proteins:
                protein.to_ttl(f)

            for isoform in isoforms:
                isoform.to_ttl(f)

            self.write_statistics(f, dataset, peptides, proteins, optimized_proteins, psms, spectra)


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

    def write_statistics(
            self, 
            f, 
            dataset: DataSet,
            peptides: list[Peptide],
            proteins: list[Protein],
            optimizaed_proteins: list[Protein],
            psms: list[Psm],
            spectra: list[Spectrum]
    ) -> None:
        f.write(f':{dataset.get_id()} \n')
        f.write('    sio:SIO_000216 [\n')
        f.write('        a jpost:NumOfPsms ;\n')
        f.write(f'        sio:SIO_000300 {len(psms)} ;\n')
        f.write('        rdfs:label "Number of PSMs"\n')
        f.write('    ] ;\n')

        f.write('    sio:SIO_000216 [\n')
        f.write('        a jpost:NumOfSpectra ;\n')
        f.write(f'        sio:SIO_000300 {len(spectra)} ;\n')
        f.write('        rdfs:label "Number of Spectra"\n')
        f.write('    ] ;\n')

        f.write('    sio:SIO_000216 [\n')
        f.write('        a jpost:NumOfPeptides ;\n')
        f.write(f'        sio:SIO_000300 {len(peptides)} ;\n')
        f.write('        rdfs:label "Number of Peptides"\n')
        f.write('    ] ;\n')

        f.write('    sio:SIO_000216 [\n')
        f.write('        a jpost:NumOfMatchedProteins ;\n')
        f.write(f'        sio:SIO_000300 {len(proteins)} ;\n')
        f.write('        rdfs:label "Number of Matched Proteins"\n')
        f.write('    ] ;\n')

        f.write('    sio:SIO_000216 [\n')
        f.write('        a jpost:NumOfProteinsWithUniquePeptide ;\n')
        f.write(f'        sio:SIO_000300 {len(optimizaed_proteins)} ;\n')
        f.write('        rdfs:label "Number of proteins with unique peptide"\n')
        f.write('    ] ;\n')

        f.write('    sio:SIO_000216 [\n')
        f.write('        a jpost:NumOfRawData ;\n')
        f.write(f'        sio:SIO_000300 {len(dataset.get_rawdata_list().get_list())} ;\n')
        f.write('        rdfs:label "Number of raw data"\n')
        f.write('    ] ;\n')


        leading_count = 0
        for protein in proteins:
            if protein.is_leading():
                leading_count += 1

        f.write('    sio:SIO_000216 [\n')
        f.write('        a jpost:NumOfLeadingProteins ;\n')
        f.write(f'        sio:SIO_000300 {leading_count} ;\n')
        f.write('        rdfs:label "Number of leading proteins"\n')
        f.write('    ] .\n\n')
