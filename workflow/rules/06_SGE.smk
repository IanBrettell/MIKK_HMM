rule sge_co_occupancy:
    input:
        data = rules.run_hmm.output,
        line_cols = rules.get_line_ranks_and_colours.output.csv,
    output:
        box_png = "book/figs/sge_F0/co-occupancy/{variables}/{interval}_{n_states}_cooc_box_all.png",
        box_pdf = "book/figs/sge_F0/co-occupancy/{variables}/{interval}_{n_states}_cooc_box_all.pdf",
        heat_png = "book/figs/sge_F0/co-occupancy/{variables}/{interval}_{n_states}_cooc_heatmap.png",
        heat_pdf = "book/figs/sge_F0/co-occupancy/{variables}/{interval}_{n_states}_cooc_heatmap.pdf",
    log:
        os.path.join(
            config["workdir"],
            "logs/sge_co_occupancy/{interval}/{variables}/{n_states}.log"
        ),
    params:
        n_states = "{n_states}",
        selected_lines = config["selected_lines"]
    resources:
        mem_mb = 30000
    container:
        config["R_4.2.0"]
    script:
        "../scripts/sge_co_occupancy.R"

rule sge_deviation:
    input:
        data = rules.run_hmm.output,
        line_cols = rules.get_line_ranks_and_colours.output.csv,
    output:
        png = "book/figs/sge_F0/ref_deviation/{variables}/{interval}_{n_states}_deviation.png",
        pdf = "book/figs/sge_F0/ref_deviation/{variables}/{interval}_{n_states}_deviation.pdf",
    log:
        os.path.join(
            config["workdir"],
            "logs/sge_deviation/{interval}/{variables}/{n_states}.log"
        ),
    params:
        n_states = "{n_states}"
    resources:
        mem_mb = 30000
    container:
        config["R_4.2.0"]
    script:
        "../scripts/sge_deviation.R"