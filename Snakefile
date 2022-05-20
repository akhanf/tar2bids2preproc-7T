from snakemake.io import glob_wildcards

configfile: 'config.yml'

globtar = glob_wildcards(config['tarfile'])

def get_tarfile(wildcards):
    subject = wildcards.subject    
    idx = globtar.subject.index(subject)

    fmt_dict = {'subject': subject}

    for tar_wc in config['tar_wildcards']:
        fmt_dict[tar_wc] = getattr(globtar,tar_wc)[idx]
    return config['tarfile'].format(**fmt_dict)

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
        participants_tsv='bids_extras/sub-{subject}/participants.tsv',
        bidsignore='bids_extras/sub-{subject}/.bidsignore',
    container: config['singularity']['tar2bids']
    shadow: 'minimal'
    shell: 
        "/opt/tar2bids/tar2bids -h {params.heuristic} "
        " -T 'sub-{{subject}}' {input} && "
        " cp -v bids/participants.tsv {output.participants_tsv} && "
        " cp -v bids/.bidsignore {output.bidsignore}"



        

