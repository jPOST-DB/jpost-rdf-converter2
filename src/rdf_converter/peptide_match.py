import os
import subprocess
from pathlib import Path
from typing import List, Optional
from .utils.logging import get_logger

logger = get_logger(__name__)

def run_peptidematch(
    jar_path: str,
    fasta_path: str,
    peptides_txt: str,
    out_tsv: str,
    java_bin: str = "java",
    extra_args: Optional[List[str]] = None,
) -> int:
    """Run the PeptideMatch command-line JAR.

    Args:
        jar_path: path to PeptideMatchCMD_1.1.jar
        fasta_path: FASTA file
        peptides_txt: A text file: one peptide per line
        out_tsv: Output TSV file path
        java_bin: Java executable
        extra_args: Additional args if needed

    Returns:
        Return code from the Java process.
    """
    Path(out_tsv).parent.mkdir(parents=True, exist_ok=True)
    cmd = [java_bin, "-jar", jar_path, "-f", fasta_path, "-p", peptides_txt, "-o", out_tsv]
    if extra_args:
        cmd.extend(extra_args)
    logger.info(f"[bold]Running PeptideMatch[/bold]: {' '.join(cmd)}")
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    logger.info(proc.stdout)
    return proc.returncode
