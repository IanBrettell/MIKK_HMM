# Calculate state frequencies for F0
rule state_freq_F2:
    input:
        data = rules.run_hmm.output,
        line_cols = rules.get_line_ranks_and_colours.output.csv
    output:
        csv_notrans = os.path.join(
            config["workdir"],
            "state_freq_F2/{variables}/{interval}/{n_states}/{dge_sge}_notrans.csv"
        ),
        csv_invnorm = os.path.join(
            config["workdir"],
            "state_freq_F2/{variables}/{interval}/{n_states}/{dge_sge}_invnorm.csv"
        ),
        hist = "book/figs/state_freq_F2/{variables}/{interval}_{n_states}_state_freq_F2_{dge_sge}.png",
    log:
        os.path.join(
            config["workdir"],
            "logs/state_freq_F0/{interval}/{variables}/{n_states}.log"
        ),
    params:
        n_states = "{n_states}",
        dge_sge = "{dge_sge}"
    resources:
        mem_mb = 80000
    container:
        config["R_4.2.0"]
    script:
        "../scripts/state_freq_F2.R"
