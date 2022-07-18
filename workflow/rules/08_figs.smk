rule karyoplots:
    input:
        data = rules.true_hmmlearn.output.csv,
        line_cols = rules.get_line_ranks_and_colours.output.csv
    output:
        karyoplot_no_missing = "book/plots/hdrr/hmmlearn_true/{bin_length}/{cov}/karyoplot_no_missing.png",
        karyoplot_with_missing = "book/plots/hdrr/hmmlearn_true/{bin_length}/{cov}/karyoplot_wi_missing.png",
    log:
        os.path.join(
            config["workdir"],
            "logs/karyoplots/{bin_length}/{cov}.log"
        ),
    params:
        cov = "{cov}",
        bin_length = "{bin_length}",
    resources:
        mem_mb = 50000
    container:
        config["R_4.1.3"]
    script:
        "../scripts/karyoplots.R"
