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
#IN_A = "/hps/nobackup/birney/users/ian/pilot/hmm_concordance_in/1/A.csv"
#IN_B = "/hps/nobackup/birney/users/ian/pilot/hmm_concordance_in/1/B.csv"
#N_STATES = int("5")
#VARIABLES = ["distance", "angle"]

## True
IN_A = snakemake.input.A
IN_B = snakemake.input.B
N_STATES = int(snakemake.params.n_states)
VARIABLES = snakemake.params.variables
OUT_FILE = snakemake.output[0]

# Set seed

np.random.seed(42)

# Read in files

dfA = pd.read_csv(IN_A)
dfB = pd.read_csv(IN_B)

# Add column indicating whether it is group A or B

dfA["group"] = "A"
dfB["group"] = "B"

# get length (number of rows) for each fish

lensA = dfA.groupby(['assay', 'date', 'time', 'quadrant', 'fish']).size().values.tolist()
lensB = dfB.groupby(['assay', 'date', 'time', 'quadrant', 'fish']).size().values.tolist()

# Select input variables

in_dfA = dfA[VARIABLES]
in_dfB = dfB[VARIABLES]

# Set up model

modelA = hmm.GaussianHMM(n_components=N_STATES, covariance_type="diag", n_iter=100)
modelB = hmm.GaussianHMM(n_components=N_STATES, covariance_type="diag", n_iter=100)

# Train

## A
modelA.fit(in_dfA, lengths = lensA)
## B
modelB.fit(in_dfB, lengths = lensB)

# Predict

statesAA = modelA.predict(in_dfA, lengths = lensA)
statesAB = modelA.predict(in_dfB, lengths = lensB)
statesBB = modelB.predict(in_dfB, lengths = lensB)
statesBA = modelB.predict(in_dfA, lengths = lensA)

# Add to dfs

dfA = dfA.assign(train_self = statesAA, train_other = statesBA)
dfB = dfB.assign(train_self = statesBB, train_other = statesAB)

# Bind into single df

df = pd.concat([dfA, dfB])

# Write to file

df.to_csv(OUT_FILE, index = False)

