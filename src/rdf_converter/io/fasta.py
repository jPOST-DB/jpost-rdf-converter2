from typing import Dict, Iterator, Tuple
from Bio import SeqIO

def read_fasta(path: str) -> Dict[str, str]:
    records = {}
    for rec in SeqIO.parse(path, "fasta"):
        records[rec.id] = str(rec.seq)
    return records

def iter_fasta(path: str) -> Iterator[Tuple[str,str]]:
    for rec in SeqIO.parse(path, "fasta"):
        yield rec.id, str(rec.seq)
