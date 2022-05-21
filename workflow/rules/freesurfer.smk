
def get_freesurfer_inputs(wildcards):
    
    checkpointrule = checkpoints.gradcorrect_extra.get(**wildcards)

    #bids exists at this point, so we can parse it 
    bids_inputs = snakebids.generate_inputs('gradcorrect',config['pybids_inputs_freesurfer'],
                participant_label=wildcards.subject,use_bids_inputs=True)

    t1_path = bids_inputs.input_path['T1w']
    subj_zip_list = snakebids.filter_list(bids_inputs.input_zip_lists['T1w'], wildcards)

    return expand(t1_path,zip,**subj_zip_list)
                            


rule freesurfer_subj:
    input:
        t1w = get_freesurfer_inputs
    params:
        in_args = lambda wildcards, input: ' '.join( f'-i {img}' for img in input.t1w ),
        subjects_dir = 'freesurfer'
    output:
        'freesurfer/sub-{subject}'
    threads: 8
    resources: 
        mem_mb=32000,
        time=1440
    shadow: 'minimal'
    container: config['singularity']['freesurfer']
    benchmark: 'benchmarks/freesurfer_sub-{subject}.tsv'
    log: 'logs/freesurfer_sub-{subject}.txt'
    shell: 
        'recon-all -threads 8 -sd {resources.tmpdir} {params.in_args} '
        '-subjid sub-{wildcards.subject} -parallel -hires -all && '
        'cp -Rv {resources.tmpdir}/sub-{wildcards.subject} {output}'
        


