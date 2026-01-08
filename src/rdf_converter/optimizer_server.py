from fastapi import FastAPI

from .optimizer_command import create_optimizer

app = FastAPI()


@app.get('/')
def help():
    return {
        'endpoints': {
            '/list': 'List of datasets',
            '/proteins': 'Get optimized proteins',
            '/peptides': 'Get optimized peptides'
        }
    }

@app.get('/list')
def list_datasets():
    '''List of datasets'''
    optimizer = create_optimizer()
    datasets = optimizer.list_datasets()
    dataset_ids = []
    for dataset in datasets:
        dataset_ids.append(dataset['dataset_id'])

    return {'datasets': dataset_ids}


@app.get('/proteins')
def get_proteins(
    dataset: str,
    min_score: float = 0.0,
    update_cache: bool = False,
    accession: bool = False
):
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

    return {'proteins': uniprots}


@app.get('/peptides')
def get_peptides(
    dataset: str,
    min_score: float = 0.0,
    update_cache: bool = False,
    accession: bool = False
):
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

            peptide_info = {'uniprot': uniprot, 'sequence': pep_seq, 'score': score}
            peptides.append(peptide_info)

    protein_map = {}
    for peptide in peptides:
        uniprot = peptide['uniprot']
        if uniprot not in protein_map:
            protein_map[uniprot] = []
        protein_map[uniprot].append({
            'sequence': peptide['sequence'],
            'score': peptide['score']
        })

    protein_list = []
    for uniprot, pep_list in protein_map.items():
        protein_list.append({
            'uniprot': uniprot,
            'peptides': pep_list
        })  

    return {'proteins': protein_list}


