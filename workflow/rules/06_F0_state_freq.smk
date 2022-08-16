# Calculate state frequencies for F0
# NOTE: this script includes writing significance tables
# to a google sheet. Therefore it needs to be run from RStudio
rule state_freq_F0:
    input:
        data = rules.run_hmm.output,
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
        sheet_id = config["aov_google_sheet_all"],
    resources:
        mem_mb = 80000
    container:
        config["R_4.2.0"]
    script:
        "../scripts/state_freq_F0.R"

# Same as above, but only with the lines selected for the F2-cross
# NOTE: this script includes writing significance tables
# to a google sheet. Therefore it needs to be run from RStudio
rule state_freq_F0_select:
    input:
        data = rules.run_hmm.output,
        line_cols = rules.get_line_ranks_and_colours.output.csv,
        selected_lines = config["selected_lines"]
    output:
        dge_hist = "book/figs/state_freq_F0_select/{variables}/{interval}_{n_states}_state_freq_F0_select_dge.png",
        sge_hist = "book/figs/state_freq_F0_select/{variables}/{interval}_{n_states}_state_freq_F0_select_sge.png"
    log:
        os.path.join(
            config["workdir"],
            "logs/state_freq_F0/{interval}/{variables}/{n_states}.log"
        ),
    params:
        n_states = "{n_states}",
        sheet_id = config["aov_google_sheet_select"],
    resources:
        mem_mb = 80000
    container:
        config["R_4.2.0"]
    script:
        "../scripts/state_freq_F0_select.R"

# Tile and time density plots for all MIKK F0
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
        sheet_id = config["aov_google_sheet_all"]
    resources:
        mem_mb = 30000
    container:
        config["R_4.2.0"]
    script:
        "../scripts/time_dependence_F0_all.R"

# Create polar plots with panels coloured by significance for lines
# selected for the F2 cross
# NOTE: this needs to be run in RStudio
rule polar_plots_with_sig:
    input:
        rules.run_hmm.output,
    output:
        polar_all_dge_of = "book/figs/polar_plots_with_sig/{interval}/{variables}/{n_states}/polar_all_dge_of.png",
        polar_all_dge_no = "book/figs/polar_plots_with_sig/{interval}/{variables}/{n_states}/polar_all_dge_no.png",
        polar_all_sge_of = "book/figs/polar_plots_with_sig/{interval}/{variables}/{n_states}/polar_all_sge_of.png",
        polar_all_sge_no = "book/figs/polar_plots_with_sig/{interval}/{variables}/{n_states}/polar_all_sge_no.png",
    log:
        os.path.join(
            config["workdir"],
            "logs/polar_plots_with_sig/{interval}/{variables}/{n_states}.log"
        ),
    params:
        interval = "{interval}",
        n_states = "{n_states}",
        sheet_id = config["aov_google_sheet_select"]
    resources:
        mem_mb = 20000
    container:
        config["R_4.2.0"]
    script:
        "../scripts/polar_plots_with_sig.R"

# Tile and time density plots for lines selected for F2 cross
rule time_dependence_F2_cross_lines:
    input:
        data = rules.split_datasets.output.F0,
        line_cols = rules.get_line_ranks_and_colours.output.csv,
    output:
        dge_of = "book/figs/time_dependence_F2_cross_lines/{variables}/select_{interval}_{n_states}_dge_of.png",
        dge_no = "book/figs/time_dependence_F2_cross_lines/{variables}/select_{interval}_{n_states}_dge_no.png",
        sge_of = "book/figs/time_dependence_F2_cross_lines/{variables}/select_{interval}_{n_states}_sge_of.png",
        sge_no = "book/figs/time_dependence_F2_cross_lines/{variables}/select_{interval}_{n_states}_sge_no.png",
    log:
        os.path.join(
            config["workdir"],
            "logs/time_dependence_F2_cross_lines/{interval}/{variables}/{n_states}.log"
        ),
    params:
        n_states = "{n_states}",
        sheet_id = lambda w: config["aov_google_sheet_select"][w.interval],
        selected_lines = config["selected_lines"]
    resources:
        mem_mb = 30000
    container:
        config["R_4.2.0"]
    script:
        "../scripts/time_dependence_F2_cross_lines.R"

