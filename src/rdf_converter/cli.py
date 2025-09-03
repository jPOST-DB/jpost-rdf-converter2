import os
from pathlib import Path
import typer
from rich import print
from dotenv import load_dotenv
from .dataset_converter import DatasetConverter
from .project_converter import ProjectConverter

app = typer.Typer(add_completion=False, help="jPOST RDF Converter (Python port)")

@app.command()
def dataset(
    tsv: str = typer.Option(..., "--tsv", help="Result TSV file"),
    fasta: str = typer.Option(..., "--fasta", help="FASTA file"),
    meta_dir: str = typer.Option(None, "--meta-dir", help="Metadata folder (optional)"),
    out: str = typer.Option(..., "--out", help="Output Turtle (TTL) path"),
    rev: str = typer.Option(..., "--rev", help="rev JPST ID"),
    branch: int = typer.Option(..., "--branch", help="Branch number"),
    pepresult: str = typer.Option(None, "--pepresult", help="PeptideMatch result TSV path"),
    peptidematch_jar: str = typer.Option(None, "--peptidematch-jar", help="Path to PeptideMatch JAR"),
    java_bin: str = typer.Option(None, "--java-bin", help="Java executable to use"),
    intermediate_dir: str = typer.Option("out", "--intermediate-dir", help="Directory for intermediate files"),
):
    """Convert dataset TSV + FASTA into TTL, with optional PeptideMatch integration."""
    load_dotenv()
    peptidematch_jar = peptidematch_jar or os.getenv("PEPTIDEMATCH_JAR")
    java_bin = java_bin or os.getenv("JAVA_BIN", "java")

    conv = DatasetConverter(tsv, fasta, meta_dir)
    conv.load()
    conv.map_columns()  # customize as needed
    conv.build_models()
    conv.write_intermediate(intermediate_dir)

    if pepresult and peptidematch_jar:
        typer.echo("Running PeptideMatch for peptides...")
        result_path = conv.run_peptidematch(peptidematch_jar, intermediate_dir, java_bin=java_bin)
        if result_path and pepresult:
            Path(pepresult).parent.mkdir(parents=True, exist_ok=True)
            # Just copy/rename produced file
            Path(result_path).replace(pepresult)
            typer.echo(f"Saved PeptideMatch results: {pepresult}")
    elif pepresult and not peptidematch_jar:
        typer.echo("[yellow]--pepresult was given but no --peptidematch-jar or PEPTIDEMATCH_JAR found; skipping[/yellow]")

    # Write minimal TTL
    dataset_id = Path(tsv).stem.upper()
    conv.to_turtle(out, dataset_id=dataset_id, rev=rev, branch=branch)
    print(f"[green]Wrote TTL:[/green] {out}")

@app.command()
def project(
    xml: str = typer.Option(..., "--xml", help="Project XML file"),
    out: str = typer.Option(..., "--out", help="Output Turtle (TTL) path"),
):
    load_dotenv()
    conv = ProjectConverter(xml)
    conv.parse()
    conv.to_turtle(out)
    print(f"[green]Wrote TTL:[/green] {out}")

if __name__ == "__main__":
    app()
