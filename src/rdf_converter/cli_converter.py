import os
from pathlib import Path
import typer
from rich import print
from dotenv import load_dotenv
from .dataset_converter import DatasetConverter

app = typer.Typer(add_completion=False, help='jPOST RDF Converter (Python port)')

@app.command()
def dataset(
    tsv: str = typer.Option(..., '--tsv', help='Result TSV file'),
    fasta: str = typer.Option(..., '--fasta', help='FASTA file'),
    meta_data: str = typer.Option(..., '--meta-data', help='Metadata'),
    out: str = typer.Option(..., '--out', help='Output Turtle (TTL) path'),
    intermediate_dir: str = typer.Option(..., '--intermediate-dir', help='Directory for intermediate files'),
    rev: str = typer.Option(..., '--rev', help='rev JPST ID'),
    pep: str = typer.Option(None, '--pep', help='PEP file (optional)'),
    branch: int = typer.Option(..., '--branch', help='Branch number'),
):
    load_dotenv()
    peptidematch_jar = os.getenv('PEPTIDEMATCH_JAR')
    java_bin = os.getenv('JAVA_BIN', 'java')

    conv = DatasetConverter(rev, branch, tsv, fasta, meta_data, pep, intermediate_dir, out, java_bin, peptidematch_jar)
    conv.convert()



if __name__ == '__main__':
    app()
