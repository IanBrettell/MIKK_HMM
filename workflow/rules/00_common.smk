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
        interval = 0.08
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 20000
    container:
        config["R_4.2.0"]
    script:
        "../scripts/get_line_ranks_and_colours.R"

