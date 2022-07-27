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
            "hmm_out_split/0.05/dist_angle/15/F0.csv",
        ),
    output:
        fig = "book/figs/line_mean_speed/line_mean_speed_0.05.png",
        csv = "config/line_colours/line_colours_0.05.csv"
    log:
        os.path.join(
            config["workdir"],
            "logs/get_line_ranks_and_colours/all.log"
        ),
    params:
        interval = 0.08,
        selected_lines = config["selected_lines"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 20000
    container:
        config["R_4.2.0"]
    script:
        "../scripts/get_line_ranks_and_colours.R"

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