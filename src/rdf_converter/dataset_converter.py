from pathlib import Path
from typing import Optional, Dict, List
import pandas as pd
from .models.core import Peptide, Protein, PSM
from .io.tsv import write_tsv
from .io.fasta import read_fasta
from .peptide_match import run_peptidematch
from .utils.logging import get_logger

logger = get_logger(__name__)

class DatasetConverter:
    """Pythonic re-implementation of Dataset conversion.

    This is a **scaffold**: adjust column mappings to your TSV format.

    """
    def __init__(self, tsv_path: str, fasta_path: str, meta_dir: Optional[str] = None):
        self.tsv_path = Path(tsv_path)
        self.fasta_path = Path(fasta_path)
        self.meta_dir = Path(meta_dir) if meta_dir else None
        self.df: Optional[pd.DataFrame] = None
        self.fasta: Dict[str, str] = {}
        self.peptides: List[Peptide] = []
        self.proteins: Dict[str, Protein] = {}
        self.psms: List[PSM] = []

    def load(self):
        logger.info(f"Reading TSV: {self.tsv_path}")
        self.df = pd.read_csv(self.tsv_path, sep='\t')
        logger.info(f"Reading FASTA: {self.fasta_path}")
        self.fasta = read_fasta(str(self.fasta_path))

    def map_columns(self, columns: Dict[str,str] = None):
        """Map your dataset's column names to standard keys.

        Provide a dict like:

            {

              'sequence': 'PeptideSequence',

              'charge': 'Charge',

              'spectrum_id': 'SpectrumID',

              'protein': 'ProteinAccession',

              'score': 'Score'

            }

        """
        if columns is None:
            columns = {
                'sequence': 'PeptideSequence',
                'charge': 'Charge',
                'spectrum_id': 'SpectrumID',
                'protein': 'ProteinAccession',
                'score': 'Score',
            }
        assert self.df is not None, "Call load() first"
        missing = [v for v in columns.values() if v not in self.df.columns]
        if missing:
            logger.warning(f"Missing columns in TSV (ok for scaffold): {missing}")
        self.df = self.df.rename(columns={v:k for k,v in columns.items() if v in self.df.columns})

    def build_models(self):
        assert self.df is not None, "Call load() first"
        for _, row in self.df.iterrows():
            seq = str(row.get('sequence', '')).strip()
            if not seq:
                continue
            pep = Peptide(sequence=seq, charge=row.get('charge'))
            prot = str(row.get('protein','')).strip() or None
            if prot:
                pep.protein_ids.append(prot)
                if prot not in self.proteins:
                    self.proteins[prot] = Protein(accession=prot)
                self.proteins[prot].peptides.append(pep)
            self.peptides.append(pep)
            self.psms.append(PSM(
                spectrum_id=str(row.get('spectrum_id','')),
                sequence=seq,
                charge=row.get('charge'),
                score=row.get('score'),
                protein=prot
            ))

    def write_intermediate(self, out_dir: str):
        out = Path(out_dir)
        write_tsv(out / "peptides.tsv", ({
            "sequence": p.sequence,
            "charge": p.charge,
            "protein_ids": ",".join(p.protein_ids),
            "modifications": p.modifications or ""
        } for p in self.peptides))

        write_tsv(out / "proteins.tsv", ({
            "accession": pr.accession,
            "n_peptides": len(pr.peptides),
            "is_leading": pr.is_leading,
            "is_anchor": pr.is_anchor
        } for pr in self.proteins.values()))

        write_tsv(out / "psms.tsv", ({
            "spectrum_id": x.spectrum_id,
            "sequence": x.sequence,
            "charge": x.charge,
            "score": x.score,
            "protein": x.protein
        } for x in self.psms))

    def run_peptidematch(self, jar_path: str, out_dir: str, java_bin: str = "java") -> Optional[str]:
        # Prepare peptides list
        out = Path(out_dir)
        peps_txt = out / "peptides.txt"
        peps_txt.parent.mkdir(parents=True, exist_ok=True)
        with open(peps_txt, "w", encoding="utf-8") as f:
            for p in sorted({p.sequence for p in self.peptides}):
                f.write(p + "\n")
        out_tsv = out / "peptidematch.tsv"
        rc = run_peptidematch(jar_path, str(self.fasta_path), str(peps_txt), str(out_tsv), java_bin=java_bin)
        if rc != 0:
            logger.warning("PeptideMatch returned non-zero exit status")
        return str(out_tsv) if out_tsv.exists() else None

    def to_turtle(self, out_ttl: str, dataset_id: str, rev: str, branch: int):
        # Minimal placeholder TTL writer. Extend as needed.
        Path(out_ttl).parent.mkdir(parents=True, exist_ok=True)
        with open(out_ttl, "w", encoding="utf-8") as w:
            w.write(f"@prefix jpost: <http://jpost.org/ontology/> .\n")
            w.write(f"@prefix : <http://jpost.org/resource/> .\n\n")
            w.write(f":{dataset_id} a jpost:Dataset ; jpost:rev \"{rev}\" ; jpost:branch {branch} .\n")
            for pr in self.proteins.values():
                w.write(f":{dataset_id} jpost:hasProtein :{pr.accession} .\n")
            for p in self.peptides:
                pid = p.sequence
                w.write(f":{dataset_id} jpost:hasPeptide :{pid} .\n")
