######################
# Libraries
######################

import pandas as pd
import numpy as np
import os
import itertools

######################
# Config file and sample sheets
######################

configfile: "config/config.yaml"

# Calculate mean speed across the course of the video for each fish
# Order lines by mean speed
# And extract colour palette for further analysis
rule get_line_ranks_and_colours:
    input:
        os.path.join(
            config["workdir"],
            "hmm_out_split/{interval}/dist_angle/15/F0.csv",
        ),
    output:
        png_all = "book/figs/line_mean_speed/line_mean_speed_{interval}_all.png",
        pdf_all = "book/figs/line_mean_speed/line_mean_speed_{interval}_all.pdf",
        png_select = "book/figs/line_mean_speed/line_mean_speed_{interval}_selected.png",
        pdf_select = "book/figs/line_mean_speed/line_mean_speed_{interval}_selected.pdf",
        csv = "config/line_colours/line_colours_{interval}.csv"
    log:
        os.path.join(
            config["workdir"],
            "logs/get_line_ranks_and_colours/{interval}.log"
        ),
    params:
        interval = "{interval}",
        selected_lines = config["selected_lines"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 20000
    container:
        config["R_4.2.0"]
    script:
        "../scripts/get_line_ranks_and_colours.R"

rule plot_line_median_and_variance:
    input:
        dat = os.path.join(
            config["workdir"],
            "hmm_out_split/{interval}/dist_angle/15/F0.csv",
        ),
        line_cols = rules.get_line_ranks_and_colours.output.csv,
    output:
        png_all = "book/figs/line_mean_speed_and_variance/line_mean_speed_variance_{interval}_all.png",
        pdf_all = "book/figs/line_mean_speed_and_variance/line_mean_speed_variance_{interval}_all.pdf",
        png_select = "book/figs/line_mean_speed_and_variance/line_mean_speed_variance_{interval}_selected.png",
        pdf_select = "book/figs/line_mean_speed_and_variance/line_mean_speed_variance_{interval}_selected.pdf",
    log:
        os.path.join(
            config["workdir"],
            "logs/plot_line_median_and_variance/{interval}.log"
        ),
    params:
        interval = 0.08,
        selected_lines = config["selected_lines"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 20000
    container:
        config["R_4.2.0"]
    script:
        "../scripts/plot_line_median_and_variance.R"


# Get the sample counts for each cross
rule count_crosses:
    input:
        config["F2_samples_file"]
    output:
        "config/F2_cross_counts.csv"
    log:
        os.path.join(
            config["workdir"],
            "logs/count_crosses/all.log"
        ),
    resources:
        mem_mb = 1000
    container:
        config["R_4.2.0"]
    script:
        "../scripts/count_crosses.R"    

######################
# Parameters
######################

PERM_SEEDS = list(range(1, config["n_permutations"][0] + 1))

# Get contigs (chromosomes 1-24)
def get_contigs(start = config["contigs"][0], end = config["contigs"][1]):
    """Get list of chromosomes."""
    end = end + 1
    return list(range(start, end))
CONTIGS = get_contigs()
