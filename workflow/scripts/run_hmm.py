import sys,os
import logging, traceback
logging.basicConfig(filename=snakemake.log[0],
                    level=logging.INFO,
                    format='%(asctime)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    )
def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logger.error(''.join(["Uncaught exception: ",
                         *traceback.format_exception(exc_type, exc_value, exc_traceback)
                         ])
                 )
# Install exception handler
sys.excepthook = handle_exception



# Import libraries
import numpy as np
import pandas as pd
from hmmlearn import hmm

# Import variables

## Debug 
#IN_FILE = "/hps/nobackup/birney/users/ian/pilot/merged/1.csv"
#N_STATES = int("5")
#VARIABLES = ["distance", "angle"]

## True
IN_FILE = snakemake.input[0]
N_STATES = int(snakemake.params.n_states)
VARIABLES = snakemake.params.variables
OUT_FILE = snakemake.output[0]

# Set seed

np.random.seed(42)

# Read in file

df = pd.read_csv(IN_FILE)

# get length (number of rows) for each sample / chr

s_chr_lens = df.groupby(['assay', 'date', 'time', 'quadrant', 'fish']).size().values.tolist()

# Select input variables

in_df = df[VARIABLES]

# Set up model

model = hmm.GaussianHMM(n_components=N_STATES, covariance_type="diag", n_iter=100)

# Train

model.fit(in_df, lengths = s_chr_lens)

# Predict

states = model.predict(in_df, lengths = s_chr_lens)

# Add to df

df = df.assign(state = states)

# Write to file

df.to_csv(OUT_FILE, index = False)

