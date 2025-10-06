# jpost_rdf_converter_cli.py
from __future__ import print_function

import os
import sys

# ---- 依存パッケージが無い環境でも動くように緩やかに import ----
try:
    from dotenv import load_dotenv  # optional
    _HAVE_DOTENV = True
except Exception:
    _HAVE_DOTENV = False

try:
    from rich import print as rprint  # optional
    def _print(msg):
        try:
            rprint(msg)
        except Exception:
            sys.stdout.write(str(msg) + "\n")
except Exception:
    def _print(msg):
        sys.stdout.write(str(msg) + "\n")

# Typer は 3.6+ 必須。無ければ argparse にフォールバック。
try:
    import typer  # optional
    _HAVE_TYPER = True
except Exception:
    _HAVE_TYPER = False

# 相対 import とスクリプト実行の両対応
try:
    # パッケージ配下（python -m package.module）のとき
    from .dataset_converter import DatasetConverter  # type: ignore
except Exception:
    # 単独スクリプト実行（python jpost_rdf_converter_cli.py）のとき
    from dataset_converter import DatasetConverter  # type: ignore


def _load_env():
    """dotenv があれば読み込む（無くても無視）。"""
    if _HAVE_DOTENV:
        try:
            load_dotenv()
        except Exception:
            pass


def run_dataset(tsv, fasta, meta_data, out_path, intermediate_dir, rev, pep, branch):
    """
    共通処理本体：CLI 実装（typer/argparse）から呼び出す。
    できるだけ古い Python でも動くよう f-string は未使用。
    """
    _load_env()

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


# ---------------- Typer（あれば使用） ----------------
if _HAVE_TYPER:
    app = typer.Typer(add_completion=False, help='jPOST RDF Converter (Python port)')

    @app.command("dataset")
    def dataset_cmd(
        tsv: str = typer.Option(..., '--tsv', help='Result TSV file'),
        fasta: str = typer.Option(..., '--fasta', help='FASTA file'),
        meta_data: str = typer.Option(..., '--meta-data', help='Metadata'),
        out: str = typer.Option(..., '--out', help='Output Turtle (TTL) path'),
        intermediate_dir: str = typer.Option(..., '--intermediate-dir', help='Directory for intermediate files'),
        rev: str = typer.Option(..., '--rev', help='rev JPST ID'),
        pep: str = typer.Option(None, '--pep', help='PEP file (optional)'),
        branch: int = typer.Option(..., '--branch', help='Branch number'),
    ):
        run_dataset(tsv, fasta, meta_data, out, intermediate_dir, rev, pep, branch)

    def main():
        app()

# ---------------- argparse フォールバック ----------------
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
        p.add_argument("--out", required=True, help="Output Turtle (TTL) path")
        p.add_argument("--intermediate-dir", required=True, help="Directory for intermediate files")
        p.add_argument("--rev", required=True, help="rev JPST ID")
        p.add_argument("--pep", required=False, default=None, help="PEP file (optional)")
        p.add_argument("--branch", required=True, type=int, help="Branch number")

        args = parser.parse_args()

        if args.command == "dataset":
            run_dataset(
                args.tsv,
                args.fasta,
                args.meta_data,
                args.out,
                args.intermediate_dir,
                args.rev,
                args.pep,
                args.branch,
            )
        else:
            # サブコマンドが無ければヘルプ表示
            parser.print_help()


if __name__ == "__main__":
    main()
