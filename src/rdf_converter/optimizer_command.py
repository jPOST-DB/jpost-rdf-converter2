import os
from pathlib import Path
import typer
from rich import print
from dotenv import load_dotenv

from rdf_converter.models import dataset
from .protein_optimizer import ProteinOptimizer
import json

app = typer.Typer(add_completion=False, help='jPOST Protein Optimizer (Python port)')


def create_optimizer():
    load_dotenv()
    api_url = os.getenv('SPARQLIST_URL', 'https://db-dev.jpostdb.org/sparqlist_pi/api')
    cache_dir = os.getenv('OPTIMIZER_CACHE', './optimizer_cache')
    optimizer = ProteinOptimizer(api_url, cache_dir)
    return optimizer


@app.command()
def list(
    output: Path = typer.Option(None, '--output', '-o', help='Output file')
):
    '''List of datasets'''
    optimizer = create_optimizer()
    datasets = optimizer.list_datasets()
    dataset_ids = []
    for dataset in datasets:
        dataset_ids.append(dataset['dataset_id'])

    if output:
        with open(output, 'w') as f:
            for dataset_id in dataset_ids:
                f.write(f'{dataset_id}\n')
    else:
        for dataset_id in dataset_ids:
            print(dataset_id)
   
   

@app.command()
def proteins(
    dataset: str = typer.Option(..., '--dataset', help='Dataset ID (comma seperated for multiple)'),
    min_score: float = typer.Option(0.0, '--min-score', help='Minimum score'),
    update_cache: bool = typer.Option(False, '--update-cache', help='Update cache'),
    output: Path = typer.Option(None, '--output', '-o', help='Output file'),
    accession: bool = typer.Option(False, '--accession', help='Show only accession (remove isoform)')
):
    '''Optimize proteins and show proteins'''
    ids = dataset.split(',')
    optimizer = create_optimizer()
    optimized_proteins = optimizer.optimize_proteins(ids, min_score, update_cache)
    uniprots = []
    for protein in optimized_proteins:
        uniprot = protein.get_uniprot()
        if accession:
            index = uniprot.find('-')
            if index > 0:
                uniprot = uniprot[:index]
        if uniprot not in uniprots:
            uniprots.append(uniprot)

    if output:
        with open(output, 'w') as f:
            for uniprot in uniprots:
                f.write(f'{uniprot}\n')
    else:
        for uniprot in uniprots:
            print(uniprot)


@app.command()
def peptides(
    dataset: str = typer.Option(..., '--dataset', help='Dataset ID (comma seperated for multiple)'),
    min_score: float = typer.Option(0.0, '--min-score', help='Minimum score'),
    update_cache: bool = typer.Option(False, '--update-cache', help='Update cache'),
    output: Path = typer.Option(None, '--output', '-o', help='Output file'),
    accession: bool = typer.Option(False, '--accession', help='Show only accession (remove isoform)')
):
    '''Optimize proteins and show proteins and peptides'''
    ids = dataset.split(',')
    optimizer = create_optimizer()
    optimized_proteins = optimizer.optimize_proteins(ids, min_score, update_cache)
    peptides = []
    for protein in optimized_proteins:
        uniprot = protein.get_uniprot()
        if accession:
            index = uniprot.find('-')
            if index > 0:
                uniprot = uniprot[:index]
        
        for match in protein.get_peptide_matches():
            peptide = match.get_peptide()
            pep_seq = peptide.get_sequence()
            score = peptide.get_score()

            peptide_info = f'{uniprot}\t{pep_seq}\t{score}'
            peptides.append(peptide_info)

    if output:
        with open(output, 'w') as f:
            for peptide in peptides:
                f.write(f'{peptide}\n')
    else:
        for peptide in peptides:
            print(peptide)





if __name__ == '__main__':
    app()
