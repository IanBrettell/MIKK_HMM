rule merge_datasets:
    input:
        expand(os.path.join(
            config["rawdir"],
            "{{interval}}/{dataset}.csv"
            ),
            dataset = config["datasets"]
        ),
    output:
        os.path.join(
            config['workdir'],
            "merged/{interval}.csv"
        ),
    log:
        os.path.join(
            config["workdir"],
            "logs/merge_datasets/{interval}.log"
        ),
    resources:
        mem_mb = 50000
    container:
        config["R_4.2.0"]
    script:
        "../scripts/merge_datasets.R"

rule run_hmm:
    input:
        rules.merge_datasets.output,
    output:
        os.path.join(
            config["workdir"],
            "hmm_out/{interval}/{variables}/{n_states}.csv"
        ),        
    log:
        os.path.join(
            config["workdir"],
            "logs/run_hmm/{interval}/{variables}/{n_states}.log"
        ),
    params:
        n_states = "{n_states}",
        variables = lambda wildcards: config["hmm_variables"][wildcards.variables]
    resources:
        # start at 5000
        mem_mb = lambda wildcards, attempt: attempt * 50000,
    container:
        config["hmmlearn"]
    script:
        "../scripts/run_hmm.py"

# Split datasets by generation (F0, F2, and Kiyosu_CC)
rule split_datasets:
    input:
        rules.run_hmm.output,
    output:
        F0 = os.path.join(
            config["workdir"],
            "hmm_out_split/{interval}/{variables}/{n_states}/F0.csv"
        ),
        F2 = os.path.join(
            config["workdir"],
            "hmm_out_split/{interval}/{variables}/{n_states}/F2.csv"
        ),
        Kiyosu_CC = os.path.join(
            config["workdir"],
            "hmm_out_split/{interval}/{variables}/{n_states}/Kiyosu_CC.csv"
        ),  
    log:
        os.path.join(
            config["workdir"],
            "logs/split_datasets/{interval}/{variables}/{n_states}.log"
        ),
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 50000,
    container:
        config["R_4.2.0"]
    script:
        "../scripts/split_datasets.R"

# We also want to test the concordance between the assigned states
# when the HMM is trained on a different dataset
# Randomise order of videos for 0.5 split into train and test datasets
rule hmm_concordance_in:
    input:
        rules.merge_datasets.output,
    output:
        A = os.path.join(
            config["workdir"],
            "hmm_concordance_in/{interval}/A.csv"
        ),
        B = os.path.join(
            config["workdir"],
            "hmm_concordance_in/{interval}/B.csv"
        ),
    log:
        os.path.join(
            config["workdir"],
            "logs/hmm_cocordance_input/{interval}.log"
        ),
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 15000
    container:
        config["rocker_tidyverse"]
    script:
        "../scripts/hmm_concordance_input.R"

# Run concordance
rule hmm_concordance_out:
    input:
        A = rules.hmm_concordance_in.output.A,
        B = rules.hmm_concordance_in.output.B,
    output:
        os.path.join(
            config["workdir"],
            "hmm_concordance_out/{interval}/{variables}/{n_states}.csv"
        ),  
    log:
        os.path.join(
            config["workdir"],
            "logs/hmm_concordance_out/{interval}/{variables}/{n_states}.log"
        ),
    params:
        n_states = "{n_states}",
        variables = lambda wildcards: config["hmm_variables"][wildcards.variables]
    resources:
        # start at 5000
        mem_mb = lambda wildcards, attempt: attempt * 15000,
    container:
        config["hmmlearn"]
    script:
        "../scripts/hmm_concordance.py"    
