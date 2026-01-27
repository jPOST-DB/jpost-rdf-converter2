# rdf-converter (Python)

A pragmatic Python3 port of the **jpost-rdf-converter** with light refactoring and a bridge to the Java **PeptideMatch** tool.

## Goals

- Provide a clean, testable Python CLI.
- Keep parity with core dataset/project conversion flows (incrementally).
- Call the bundled `PeptideMatchCMD_1.0.jar` from Python.
- Emit **intermediate CSV/TSV** artifacts for traceability (peptides, proteins, PSMs).

## Quick start

```bash
# (Optionally) create a venv
python -m venv .venv && source .venv/bin/activate  # Windows: .venv\Scripts\activate

# Upgrade pip
python -m pip install --upgrade pip

# Install
pip install -r requirements.txt
pip install -e .

# Show help
rdf-convert --help

# Example (RDF化: dataset)
rdf-convert dataset \
  --tsv example/results.tsv \
  --fasta example/proteins.fasta \
  --meta_data example/meta.xml \
  --out out/example.ttl \
  --intermediate_dir out \
  --rev JPST000000 \
  --branch 1 \
  --pep pep.txt

# Example (RDF化: project)
rdf-convert project --meta_data example/project.xml --out out/project.ttl --rev JPST000000

# Example (最適化: list)
protein-optimize list

# Example (最適化: server)
protein-optimize server --port 8081

# Example (最適化: proteins)
protein-optimize proteins --dataset DS0000_1 --min-score --output out/proteins.txt

# Example (最適化: peptides)
protein-optimize peptides --dataset DS0000_1 --min-score --output out/peptides.txt

```

## Intermediate artifacts

- `out/peptides.tsv`
- `out/proteins.tsv`
- `out/psms.tsv`
- `out/peptidematch.tsv` (when running PeptideMatch)

## Configuration

You can place a `.env` file at the project root with:

```
PEPTIDEMATCH_JAR=./lib/PeptideMatchCMD_1.1.jar
JAVA_BIN=java
```

## License

This scaffold is provided for internal porting. Verify third-party licenses (e.g., PeptideMatch).

