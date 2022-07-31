####################
# Codon cluster
####################

ssh codon
module load singularity-3.7.0-gcc-9.3.0-dp5ffrp
bsub -Is bash
# If needing to copy videos from FTP (rule copy_videos),
# Need to use the datamover queue so that it can see the FTP drive:
# bsub -M 20000 -q datamover -Is bash
cd /hps/software/users/birney/ian/repos/MIKK_HMM
conda activate snakemake_7.8.2
# Regular
snakemake \
  --jobs 5000 \
  --latency-wait 300 \
  --cluster-config config/cluster.yaml \
  --cluster 'bsub -g /snakemake_bgenie -J {cluster.name} -q {cluster.queue} -n {cluster.n} -M {cluster.memory} -o {cluster.outfile}' \
  --keep-going \
  --rerun-incomplete \
  --use-conda \
  --use-singularity \
  --rerun-triggers mtime \
  --restart-times 0 \
  -s workflow/Snakefile \
  -p

# R container build
bsub -M 20000 -Is bash
module load singularity-3.7.0-gcc-9.3.0-dp5ffrp
cd /hps/software/users/birney/ian/repos/MIKK_HMM
RCONT=/hps/nobackup/birney/users/ian/containers/MIKK_HMM/R_4.2.0.sif
singularity build --remote \
    $RCONT \
    workflow/envs/R_4.2.0/R_4.2.0.def

# RStudio container run
ssh proxy-codon
bsub -M 50000 -q short -Is bash
module load singularity-3.7.0-gcc-9.3.0-dp5ffrp
cd /hps/software/users/birney/ian/repos/MIKK_HMM
CONT=/hps/nobackup/birney/users/ian/containers/MIKK_HMM/R_4.2.0.sif
singularity shell --bind /hps/nobackup/birney/users/ian/R_tmp/R_4.2.0/rstudio_db:/var/lib/rstudio-server \
                  --bind /hps/nobackup/birney/users/ian/R_tmp/R_4.2.0/tmp:/tmp \
                  --bind /hps/nobackup/birney/users/ian/R_tmp/R_4.2.0/run:/run \
                  $CONT
rstudio-server kill-all
rserver \
    --rsession-config-file /hps/software/users/birney/ian/repos/MIKK_HMM/workflow/envs/R_4.2.0/rsession.conf \
    --server-user brettell

ssh -L 8787:hl-codon-37-04:8787 proxy-codon

# RStudio for KaryoploteR
ssh proxy-codon
bsub -M 50000 -q short -Is bash
module load singularity-3.7.0-gcc-9.3.0-dp5ffrp
cd /hps/software/users/birney/ian/repos/MIKK_HMM
CONT=/hps/nobackup/birney/users/ian/containers/MIKK_HMM/R_4.1.3.sif
singularity shell --bind /hps/nobackup/birney/users/ian/R_tmp/R_4.1.3/rstudio_db:/var/lib/rstudio-server \
                  --bind /hps/nobackup/birney/users/ian/R_tmp/R_4.1.3/tmp:/tmp \
                  --bind /hps/nobackup/birney/users/ian/R_tmp/R_4.1.3/run:/run \
                  $CONT
rstudio-server kill-all
rserver \
    --rsession-config-file /hps/software/users/birney/ian/repos/MIKK_HMM/workflow/envs/R_4.1.3/rsession.conf \
    --server-user brettell

ssh -L 8787:hl-codon-37-04:8787 proxy-codon


# Faspex
## Download aspera
mamba install -c hcc aspera-cli
# Get URL
fURL=`aspera faspex list --user="felix.loosli@kit.edu" \
                         --host="faspex.embl.de" \
                         --password='MIKKF2#' \
                         --insecure | grep "ˆFaspex URL: " | awk '{ print $NF }'`
