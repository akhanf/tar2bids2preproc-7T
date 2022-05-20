from snakemake.io import glob_wildcards

configfile: 'config.yml'

globtar = glob_wildcards(config['tarfile'])

subjects = globtar.subject

def get_tarfile(wildcards):
    subject = wildcards.subject    
    idx = globtar.subject.index(subject)

    fmt_dict = {'subject': subject}

    for tar_wc in config['tar_wildcards']:
        fmt_dict[tar_wc] = getattr(globtar,tar_wc)[idx]
    return config['tarfile'].format(**fmt_dict)

rule all_bids:
    input: 
        'bids/code/validator_output.txt',


rule link_tarfile:
    input:
        get_tarfile
    output:
        'tar/sub-{subject}.tar'
    shell:
        'ln -srv {input} {output}'

rule tar2bids:
    input: rules.link_tarfile.output
    params:
        heuristic='cfmm_bold_rest.py'
    output: 
        subject_dir=directory('bids/sub-{subject}'),
        participants_tsv=temp('bids/sub-{subject}_participants.tsv'),
        bidsignore=temp('bids/sub-{subject}_bidsignore'),
        dd=temp('bids/sub-{subject}_dataset_description.json'),
    container: config['singularity']['tar2bids']
    shadow: 'minimal'
    shell: 
        "/opt/tar2bids/tar2bids -h {params.heuristic} "
        " -T 'sub-{{subject}}' {input} && "
        " cp -v bids/participants.tsv {output.participants_tsv} && "
        " cp -v bids/.bidsignore {output.bidsignore} && "
        " cp -v bids/dataset_description.json {output.dd}"



rule merge_bidsignore:
    """just gets the first bidsignore, since safe to assume all will be the same"""
    input: 
        bidsignore=expand(rules.tar2bids.output.bidsignore,subject=subjects),
    output:
        bidsignore='bids/.bidsignore'
    shell:
        'cp {input[0]} {output} '

rule merge_datasetdescription:
    """just gets the first dataset_description, since safe to assume all will be the same
        TODO: create this from config.yml"""

    input: 
        dataset_description=expand(rules.tar2bids.output.dd,subject=subjects),
    output:
        dataset_description='bids/dataset_description.json'
    shell:
        'cp {input[0]} {output} '


rule merge_participants_tsv:
    input: 
        participants_tsv=expand(rules.tar2bids.output.participants_tsv,subject=subjects),
    output:
        participants_tsv='bids/participants.tsv'
    shell:
        'echo participant_id > {output} && '
        'grep -h sub {input} | uniq | sort >> {output}'
    
rule validator:
    input: 
        'bids/participants.tsv',
        'bids/.bidsignore',
        'bids/dataset_description.json',
    output: 
        'bids/code/validator_output.txt'
    container: config['singularity']['tar2bids']
    shell: 'bids-validator bids | tee {output}'



