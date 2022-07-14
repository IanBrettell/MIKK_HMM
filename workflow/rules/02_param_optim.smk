# Recode the states for group B so that they're effectively the same state as group A
rule recode_concordance:
    input:
        rules.hmm_concordance_out.output,
    output:
        rds = os.path.join(
            config["workdir"],
            "hmm_concordance_recode/{interval}/{variables}/{n_states}.rds"
        ),
        confmat_png = "book/figs/hmm_concordance/{interval}/{variables}/{n_states}/conf_mat.png",
        confmat_pdf = "book/figs/hmm_concordance/{interval}/{variables}/{n_states}/conf_mat.pdf",
        polar_png = "book/figs/hmm_concordance/{interval}/{variables}/{n_states}/eye_facet.png",
        polar_pdf = "book/figs/hmm_concordance/{interval}/{variables}/{n_states}/eye_facet.pdf"
    log:
        os.path.join(
            config["workdir"],
            "logs/recode_concordance/{interval}/{variables}/{n_states}.log"
        ),
    params:
        interval = "{interval}",
        n_states = "{n_states}",
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 20000,
    container:
        config["R"]
    script:
        "../scripts/recode_concordance.R"

rule run_kruskal_wallis:
    input:
        expand(rules.split_datasets.output.F0,
                interval = config["seconds_interval"],
                variables = config["hmm_variables"],
                n_states = config["n_states"]               
        ),
    output:
        os.path.join(
            config["workdir"],
            "kruskal_wallis/out.rds"
        ),   
    log:
        os.path.join(
            config["workdir"],
            "logs/run_kruskal_wallis/all.log"
        ),
    resources:
        # start at 50000
        mem_mb = lambda wildcards, attempt: attempt * 80000
    container:
        config["R"]
    script:
        "../scripts/run_kruskal_wallis.R"

rule compare_params:
    input:
        conc = expand(rules.recode_concordance.output.rds,
            interval = config["seconds_interval"],
            variables = config["hmm_variables"],
            n_states = config["n_states"]
        ),
        kw = rules.run_kruskal_wallis.output,
    output:
        png = "book/figs/compare_params/compare_params.png",
        pdf = "book/figs/compare_params/compare_params.pdf"
    log:
        os.path.join(
            config["workdir"],
            "logs/compare_params/all.log"
        ),
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 20000,
    container:
        config["R_4.2.0"]
    script:
        "../scripts/compare_params.R"

