from dataclasses import dataclass, field
from typing import List, Optional, Dict

@dataclass
class Peptide:
    sequence: str
    charge: Optional[int] = None
    modifications: Optional[str] = None
    protein_ids: List[str] = field(default_factory=list)

@dataclass
class Protein:
    accession: str
    description: Optional[str] = None
    length: Optional[int] = None
    gene: Optional[str] = None
    is_leading: bool = False
    is_anchor: bool = False
    peptides: List[Peptide] = field(default_factory=list)
    meta: Dict[str, str] = field(default_factory=dict)

@dataclass
class PSM:
    spectrum_id: str
    sequence: str
    charge: Optional[int] = None
    score: Optional[float] = None
    protein: Optional[str] = None
