
rule fmriprep_subj:
    input:
        gradcorrect=get_gradcorrect,
        freesurfer='freesurfer/sub-{subject}'
    params:
        fs_subjects_dir='freesurfer',
        fs_license_file=config['fs_license_file']
    output:
        directory('fmriprep/sub-{subject}'),
        dd='fmriprep-extra/sub-{subject}_dataset_description.json'
    container: config['singularity']['fmriprep']
    shadow: 'minimal'
    benchmark: 'benchmarks/fmriprep_sub-{subject}.tsv'
    log: 'logs/fmriprep_sub-{subject}.txt'
    threads: 8
    resources: 
        mem_mb=32000,
        time=360
    shell: 
        'fmriprep {input.gradcorrect} fmriprep participant --participant_label {wildcards.subject} '
        '--nthreads {threads} --n_cpus {threads} --mem_mb {resources.mem_mb} --omp-nthreads {threads} '
        '--fs-license-file {params.fs_license_file} --fs-subjects-dir {params.fs_subjects_dir} --cifti-output '
        '--notrack --output-layout bids --use-aroma '
        '--output-spaces T1w MNI152NLin2009cAsym MNI152NLin6Asym && ' 
        'cp fmriprep/dataset_description.json {output.dd}'


checkpoint fmriprep_extra:
    input:
        dd=expand(rules.fmriprep_subj.output.dd,subject=subjects)
    output:
        dd='fmriprep/dataset_description.json'
    shell:
        'cp {input[0]} {output}'

 
