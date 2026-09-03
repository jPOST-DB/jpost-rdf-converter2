import os
from pathlib import Path
import typer
from dotenv import load_dotenv
from .dataset_converter import DatasetConverter
from .project_converter import ProjectConverter
from .models.project import Project

app = typer.Typer(add_completion=False, help='jPOST RDF Converter (Python port)')


def _resolve_rev(rev: str | None, meta_data: str) -> str:
    if rev:
        return rev
    project = Project.read_project(meta_data)
    resolved = project.get_id()
    if not resolved:
        raise typer.BadParameter('Could not detect JPST ID from metadata. Please specify --rev explicitly.')
    typer.echo(f'Auto-detected JPST ID: {resolved}')
    return resolved


@app.command()
def dataset(
    tsv: str = typer.Option(..., '--tsv', help='Result TSV file'),
    fasta: str = typer.Option(..., '--fasta', help='FASTA file'),
    meta_data: str = typer.Option(..., '--meta-data', help='Metadata'),
    output: str = typer.Option(..., '--output', '-o', help='Output Turtle (TTL) path'),
    intermediate_dir: str = typer.Option(..., '--intermediate-dir', help='Directory for intermediate files'),
    rev: str = typer.Option(None, '--rev', help='rev JPST ID (auto-detected from metadata if omitted)'),
    pep: str = typer.Option(None, '--pep', help='PEP file (optional)'),
    branch: int = typer.Option(..., '--branch', help='Branch number'),
):
    load_dotenv()
    peptidematch_jar = os.getenv('PEPTIDEMATCH_JAR')
    java_bin = os.getenv('JAVA_BIN', 'java')

    rev = _resolve_rev(rev, meta_data)
    conv = DatasetConverter(rev, branch, tsv, fasta, meta_data, pep, intermediate_dir, output, java_bin, peptidematch_jar)
    conv.convert()

@app.command()
def project(
    meta_data: str = typer.Option(..., '--meta-data', help='Metadata'),
    rev: str = typer.Option(None, '--rev', help='rev JPST ID (auto-detected from metadata if omitted)'),
    output: str = typer.Option(..., '--output', '-o', help='Output Turtle (TTL) path'),
):
    load_dotenv()
    rev = _resolve_rev(rev, meta_data)
    conv = ProjectConverter(rev, meta_data, output)
    conv.convert()


if __name__ == '__main__':
    app()
