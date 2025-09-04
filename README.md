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

# Install
pip install -r requirements.txt
pip install -e .

# Show help
rdf-convert -help

# Example (dataset) — flags are compatible conceptually, not 1:1 yet
rdf-convert dataset \
  --tsv example/results.tsv \
  --fasta example/proteins.fasta \
  --meta-dir example/meta \
  --out out/example.ttl \
  --rev JPST000000 \
  --branch 1 \
  --pepresult out/peptidematch.tsv

# Example (project)
jpconvert project --xml example/project.xml --out out/project.ttl
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

