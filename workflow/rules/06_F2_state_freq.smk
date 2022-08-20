# Calculate state frequencies for F0
rule state_freq_F2:
    input:
        data = rules.run_hmm.output,
    output:
        csv_notrans = expand(os.path.join(
            config["workdir"],
            "state_freq_F2/{{variables}}/{{interval}}/{{n_states}}/{{dge_sge}}/notrans/{state}.csv"
            ),
                state = list(range(1,16))
        ),
        csv_invnorm = expand(os.path.join(
            config["workdir"],
            "state_freq_F2/{{variables}}/{{interval}}/{{n_states}}/{{dge_sge}}/invnorm/{state}.csv"
            ),
                state = list(range(1,16))
        ),
        hist_png = "book/figs/state_freq_F2/{variables}/{interval}_{n_states}_state_freq_F2_{dge_sge}.png",
        hist_pdf = "book/figs/state_freq_F2/{variables}/{interval}_{n_states}_state_freq_F2_{dge_sge}.pdf",
    log:
        os.path.join(
            config["workdir"],
            "logs/state_freq_F2/{interval}/{variables}/{n_states}/{dge_sge}.log"
        ),
    params:
        dge_sge = "{dge_sge}"
    resources:
        mem_mb = 80000
    container:
        config["R_4.2.0"]
    script:
        "../scripts/state_freq_F2.R"

# Calculate each F2's mean speed
rule mean_speed_F2:
    input:
        data = expand(os.path.join(
            config["workdir"],
            "hmm_out_split/{{interval}}/{variables}/{n_states}/F2.csv"),
                variables = "dist_angle",
                n_states = 15,
        ),
    output:
        csv_notrans = os.path.join(
            config["workdir"],
            "mean_speed_F2/{interval}/{dge_sge}/notrans.csv"
            ),
        csv_invnorm = os.path.join(
            config["workdir"],
            "mean_speed_F2/{interval}/{dge_sge}/invnorm.csv"
            ),
        hist = "book/figs/state_freq_F2/{interval}/{dge_sge}.png",
    log:
        os.path.join(
            config["workdir"],
            "logs/state_freq_F2/{interval}/{dge_sge}.log"
        ),
    params:
        dge_sge = "{dge_sge}"
    resources:
        mem_mb = 80000
    container:
        config["R_4.2.0"]
    script:
        "../scripts/mean_speed_F2.R"

