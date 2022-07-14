rule get_line_ranks_and_colours:
    input:
        rules.merge_datasets.output,
    output:
        fig = "book/figs/line_mean_speed/line_mean_speed_{interval}.png",
        csv = "config/line_colours/line_colours_{interval}.csv"
    log:
        os.path.join(
            config["working_dir"],
            "logs/get_line_ranks_and_colours/{interval}.log"
        ),
    params:
        interval = "{interval}"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 15000
    container:
        config["R_4.2.0"]
    script:
        "../scripts/get_line_ranks_and_colours.R"
