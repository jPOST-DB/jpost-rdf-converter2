from __future__ import print_function

import os
import sys

try:
    from dotenv import load_dotenv  # optional
    _HAVE_DOTENV = True
except Exception:
    _HAVE_DOTENV = False

try:
    import typer  # optional
    _HAVE_TYPER = True
except Exception:
    _HAVE_TYPER = False

try:
    from .dataset_converter import DatasetConverter  # type: ignore
except Exception:
    from dataset_converter import DatasetConverter  # type: ignore


def _load_env():
    if _HAVE_DOTENV:
        try:
            load_dotenv()
        except Exception:
            pass


def _resolve_rev(rev, meta_data):
    if rev:
        return rev
    try:
        import xml.etree.ElementTree as ET
        tree = ET.parse(meta_data)
        root = tree.getroot()
        node = root.find('Project')
        resolved = node.get('id') if node is not None else None
        if resolved:
            _print("Auto-detected JPST ID: {0}".format(resolved))
            return resolved
    except Exception:
        pass
    raise SystemExit('ERROR: Could not detect JPST ID from metadata. Please specify --rev explicitly.')


def _print(msg):
    print(msg)


def run_dataset(tsv, fasta, meta_data, out_path, intermediate_dir, rev, pep, branch):
    _load_env()

    rev = _resolve_rev(rev, meta_data)

    peptidematch_jar = os.getenv('PEPTIDEMATCH_JAR')
    java_bin = os.getenv('JAVA_BIN', 'java')

    _print("[jPOST RDF Converter] starting convert...")
    _print("  rev={0}, branch={1}".format(rev, branch))
    _print("  tsv={0}".format(tsv))
    _print("  fasta={0}".format(fasta))
    _print("  meta_data={0}".format(meta_data))
    if pep:
        _print("  pep={0}".format(pep))
    _print("  intermediate_dir={0}".format(intermediate_dir))
    _print("  out={0}".format(out_path))
    if peptidematch_jar:
        _print("  PEPTIDEMATCH_JAR={0}".format(peptidematch_jar))
    _print("  JAVA_BIN={0}".format(java_bin))

    conv = DatasetConverter(
        rev,
        branch,
        tsv,
        fasta,
        meta_data,
        pep,
        intermediate_dir,
        out_path,
        java_bin,
        peptidematch_jar
    )
    conv.convert()
    _print("[jPOST RDF Converter] done.")


if _HAVE_TYPER:
    app = typer.Typer(add_completion=False, help='jPOST RDF Converter (Python port)')

    @app.command("dataset")
    def dataset_cmd(
        tsv: str = typer.Option(..., '--tsv', help='Result TSV file'),
        fasta: str = typer.Option(..., '--fasta', help='FASTA file'),
        meta_data: str = typer.Option(..., '--meta-data', help='Metadata'),
        output: str = typer.Option(..., '--output', '-o', help='Output Turtle (TTL) path'),
        intermediate_dir: str = typer.Option(..., '--intermediate-dir', help='Directory for intermediate files'),
        rev: str = typer.Option(None, '--rev', help='rev JPST ID (auto-detected from metadata if omitted)'),
        pep: str = typer.Option(None, '--pep', help='PEP file (optional)'),
        branch: int = typer.Option(..., '--branch', help='Branch number'),
    ):
        run_dataset(tsv, fasta, meta_data, output, intermediate_dir, rev, pep, branch)

    def main():
        app()
else:
    import argparse

    def main():
        parser = argparse.ArgumentParser(
            prog="jpost-rdf-converter",
            description="jPOST RDF Converter (Python port)"
        )
        sub = parser.add_subparsers(dest="command")

        p = sub.add_parser("dataset", help="Convert dataset to TTL")
        p.add_argument("--tsv", required=True, help="Result TSV file")
        p.add_argument("--fasta", required=True, help="FASTA file")
        p.add_argument("--meta-data", required=True, help="Metadata")
        p.add_argument("--output", "-o", required=True, help="Output Turtle (TTL) path")
        p.add_argument("--intermediate-dir", required=True, help="Directory for intermediate files")
        p.add_argument("--rev", required=False, default=None, help="rev JPST ID (auto-detected from metadata if omitted)")
        p.add_argument("--pep", required=False, default=None, help="PEP file (optional)")
        p.add_argument("--branch", required=True, type=int, help="Branch number")

        args = parser.parse_args()

        if args.command == "dataset":
            run_dataset(
                args.tsv,
                args.fasta,
                args.meta_data,
                args.output,
                args.intermediate_dir,
                args.rev,
                args.pep,
                args.branch,
            )
        else:
            parser.print_help()


if __name__ == "__main__":
    main()
