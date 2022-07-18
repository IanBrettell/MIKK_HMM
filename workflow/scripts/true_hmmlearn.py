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
import pickle

# Import variables

## Debug
IN = "/hps/nobackup/birney/users/ian/MIKK_HMM/hmm_in/F2/hdrr/20000.csv"

## True
IN = snakemake.input[0]
COV = float(snakemake.params.cov)
#LOW_COV_SAMPLES = snakemake.params.low_cov_samples
#LOW_COV_SAMPLES = [int(i) for i in LOW_COV_SAMPLES] # convert to integers
OUT_CSV = snakemake.output.csv
OUT_PCK = snakemake.output.pck

######################
# Setup
######################

# Read in file

df = pd.read_csv(IN, comment = "#")

# get length (number of rows) for each sample / chr

s_chr_lens = df.groupby(['SAMPLE', 'CHROM']).size().values.tolist()


#######################
## Filter DF for first five samples (hash out when running on full sample)
#######################
#
#s_keep = df['SAMPLE'].unique()[0:10].tolist()
#df = df[df['SAMPLE'].isin(s_keep)]
#
#s_chr_lens = df.groupby(['SAMPLE', 'CHROM']).size().values.tolist()

# Create HMM input

hmm_in = df.loc[:, ['PROP_MAT']] # Needs to remain as a data frame

######################
# HMM mods
######################

# Set initial means so that the states are in the right order (0,1,2 for Cab, Het, Kaga)
means = np.array([[0],
                 [0.5],
                 [1]])

# Params: GaussianHMM(params='stmc', init_params='stmc')
# From the docs:
# The parameters that get updated during (params) or initialized before (init_params) the training. 
# Can contain any combination of ‘s’ for startprob, ‘t’ for transmat, and other characters for subclass-specific emission parameters. 
# Defaults to all parameters.

# So if you want to initialize the means, set `init_params = "stc"`

base_model = hmm.GaussianHMM(n_components=3,
                             covariance_type="diag",
                             n_iter=100,
                             algorithm = 'viterbi',
                             init_params = 'stc' # 's' (startprob), 't' (transmat), 'm' (means), or 'c' (covars)
                             )
                             
transmat = np.array([[0.999, 0.00066, 0.00033],
                     [0.0005, 0.999, 0.0005],
                     [0.00033, 0.00066, 0.999]
                     ])


model = base_model
model.means_ = means
model.fit(hmm_in, lengths = s_chr_lens)
model.transmat_ = transmat
# Note: in the models above, these tend to be 
# [0.001, 0.09, 0.001]
model.covars_ = np.array([
                   [COV],   # HOM_PAT true 
                   [COV],   # HET true
                   [COV],   # HOM_MAT true 
                  ])
states = model.predict(hmm_in, lengths = s_chr_lens)


######################
# Add to df
######################

df = df.assign(STATE = states)

######################
# Write to file
######################

## Model
pickle.dump(model, open(OUT_PCK, 'wb'))
## csv
df.to_csv(OUT_CSV, index = False)

