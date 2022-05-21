localrules: link_tarfile,merge_bidsignore,merge_dataset_description,merge_participants_tsv,validator

def get_tarfile(wildcards):
    subject = wildcards.subject    

    idx = globtar.subject.index(subject)
    fmt_dict = dict()
    for tar_wc in get_wildcard_names(config['tarfile']):
        fmt_dict[tar_wc] = getattr(globtar,tar_wc)[idx]

    return strip_wildcard_constraints(config['tarfile']).format(**fmt_dict)


rule link_tarfile:
    input:
        get_tarfile
    output:
        'tar/sub-{subject}.tar'
    shell:
        'ln -srv {input} {output}'

rule tar2bids:
    input: 
        tarfile=rules.link_tarfile.output,
        heuristic=config['heuristic']
    params:
    output: 
        subject_dir=directory('bids/sub-{subject}'),
        participants_tsv='bids-extra/sub-{subject}_participants.tsv',
        bidsignore='bids-extra/sub-{subject}_bidsignore',
        dd='bids-extra/sub-{subject}_dataset_description.json',
    container: config['singularity']['tar2bids']
    log: 'logs/tar2bids_sub-{subject}.txt'
    benchmark: 'benchmarks/tar2bids_sub-{subject}.tsv'
    shadow: 'minimal'
    threads: 8
    resources: 
        mem_mb=32000,
        time=60
    shell: 
        "/opt/tar2bids/tar2bids -h {input.heuristic} "
        " -T 'sub-{{subject}}' {input.tarfile} && "
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

rule merge_participants_tsv:
    input: 
        participants_tsv=expand(rules.tar2bids.output.participants_tsv,subject=subjects),
    output:
        participants_tsv='bids/participants.tsv'
    shell:
        'echo participant_id > {output} && '
        'grep -h sub {input} | sort | uniq >> {output}'
 
rule merge_dataset_description:
    """just gets the first dataset_description, since safe to assume all will be the same
        TODO: create this from config.yml"""

    input: 
        bidsignore='bids/.bidsignore',
        participants_tsv='bids/participants.tsv',
        dataset_description=expand(rules.tar2bids.output.dd,subject=subjects),
    output:
        dataset_description='bids/dataset_description.json'
    shell:
        'cp {input[0]} {output} '

   
checkpoint validator:
    input: 
        'bids/dataset_description.json',
    output: 
        'bids/code/validator_output.txt'
    container: config['singularity']['tar2bids']
    shell: 'bids-validator bids | tee {output}'

def get_bids(wildcards):
    checkpointrule = checkpoints.validator.get(**wildcards)
    return 'bids'

 
