# Calculate state frequencies for F0
rule state_freq_F0:
    input:
        data = rules.split_datasets.output.F0,
        line_cols = rules.get_line_ranks_and_colours.output.csv
    output:
        dge_hist = "book/figs/state_freq_F0/{variables}/{interval}_{n_states}_state_freq_F0_dge.png",
        sge_hist = "book/figs/state_freq_F0/{variables}/{interval}_{n_states}_state_freq_F0_sge.png"
    log:
        os.path.join(
            config["workdir"],
            "logs/state_freq_F0/{interval}/{variables}/{n_states}.log"
        ),
    params:
        n_states = "{n_states}",
        sheet_id = config["aov_google_sheet"]
    resources:
        mem_mb = 80000
    container:
        config["R_4.2.0"]
    script:
        "../scripts/state_freq_F0.R"

# Tile, time density, and spatial density plots for all MIKK F0
rule time_dependence_F0_all:
    input:
        data = rules.split_datasets.output.F0,
        line_cols = rules.get_line_ranks_and_colours.output.csv
    output:
        tile_dge = "book/figs/time_dependence_F0/{variables}/{interval}_{n_states}_tile_dge.png",
        tile_sge = "book/figs/time_dependence_F0/{variables}/{interval}_{n_states}_tile_sge.png",
        dens_dge = "book/figs/time_dependence_F0/{variables}/{interval}_{n_states}_dens_dge.png",
        dens_sge = "book/figs/time_dependence_F0/{variables}/{interval}_{n_states}_dens_sge.png",
    log:
        os.path.join(
            config["workdir"],
            "logs/time_dependence_F0_all/{interval}/{variables}/{n_states}.log"
        ),
    params:
        n_states = "{n_states}",
        sheet_id = config["aov_google_sheet"]
    resources:
        mem_mb = 10000
    container:
        config["R_4.2.0"]
    script:
        "../scripts/time_dependence_F0_all.R"
