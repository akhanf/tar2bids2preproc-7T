rule gradcorrect_subj:
    input:
        get_bids,
        coeff=config['grad_coeff_file']
    output:
        directory('gradcorrect/sub-{subject}')
    shadow:
        'minimal'
    container: config['singularity']['gradcorrect']
    log: 'logs/gradcorrect_sub-{subject}.txt'
    benchmark: 'benchmarks/gradcorrect_sub-{subject}.tsv'
    threads: 4
    resources: 
        mem_mb=16000,
        time=180
    shell:
        '/gradcorrect/run.sh bids gradcorrect participant --participant_label {wildcards.subject} --grad_coeff_file {input.coeff} &> {log}'
        

checkpoint gradcorrect_extra:
    input:
        expand('gradcorrect/sub-{subject}',subject=subjects),
        bidsfiles=multiext('bids/','participants.tsv','.bidsignore','dataset_description.json')
    output:
        bidsfiles=multiext('gradcorrect/','participants.tsv','.bidsignore','dataset_description.json')
    run:
        for i,o in zip(input.bidsfiles,output.bidsfiles):
            shell('cp {i} {o}')
        
def get_gradcorrect(wildcards):
    checkpointrule = checkpoints.gradcorrect_extra.get(**wildcards)
    return 'gradcorrect'



