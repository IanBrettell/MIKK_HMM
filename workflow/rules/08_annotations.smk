## pull all significant SNPs into a file for annotation
## also add REF and ALT alleles
#rule pull_significant_snps_sf:
#    input:
#        # MLMA output
#        res = expand(os.path.join(
#            config["workdir"],
#            "gcta/mlma_loco_state_freq_consol/true/hdrr/{variables}/{{interval}}/{n_states}/{bin_length}/{cov}/{dge_sge}/{transformation}/{assay}/{state}/{covars}.mlma"
#            ),
#                variables = "dist_angle",
#                n_states = 15,
#                bin_length = 5000,
#                cov = 0.8,
#                dge_sge = ["dge", "sge"],
#                transformation = "invnorm",
#                assay = ["open_field", "novel_object"],
#                state = list(range(1,16)),
#                covars = "time-quadrant"
#        ),
#        # Minimum p-value per permutation
#        min_p = expand(os.path.join(
#            config["workdir"],
#            "gcta/mlma_loco/min_p/hdrr/{variables}/{{interval}}/{n_states}/{bin_length}/{cov}/{dge_sge}/{transformation}/{assay}/{state}/{covars}.csv"
#            ),
#                variables = "dist_angle",
#                n_states = 15,
#                bin_length = 5000,
#                cov = 0.8,
#                dge_sge = ["dge", "sge"],
#                transformation = "invnorm",
#                assay = ["open_field", "novel_object"],
#                state = list(range(1,16)),
#                covars = "time-quadrant"                
#        ),
#        # All bi-allelic alleles (REF/ALT) in the F0 .vcf
#        alleles = rules.extract_all_snps.output,
#    output:
#        # Significant SNPs
#        sigs = os.path.join(
#            config["workdir"],
#            "sig_snps/hdrr/dist_angle/{interval}/15/5000/0.8/invnorm/sigs.csv"
#        ),
#    log:
#        os.path.join(
#            config["workdir"],
#            "logs/pull_significant_snps_sf/hdrr/dist_angle/{interval}/15/5000/0.8/invnorm.log"
#        ),
#    resources:
#        mem_mb = 80000
#    container:
#        config["R_4.2.0"]
#    script:
#        "../scripts/pull_significant_snps_sf.R"
#
#rule prepare_vep_input:
#    input:
#        rules.pull_significant_snps_sf.output.sigs,
#    output:
#        out = "results/annotations/hdrr/dist_angle/{interval}/15/5000/0.8/invnorm/{interval}_vep_in.txt",
#    log:
#        os.path.join(
#            config["workdir"],
#            "logs/prepare_vep_input/hdrr/dist_angle/{interval}/15/5000/0.8/invnorm.log"
#        ),
#    resources:
#        mem_mb = 500
#    container:
#        config["tidyverse_4.1.3"]
#    script:
#        "../scripts/prepare_vep_input.R"
#
#rule run_vep_invnorm:
#    input:
#        rules.prepare_vep_input.output.out,
#    output:
#        out = "results/annotations/hdrr/dist_angle/{interval}/15/5000/0.8/invnorm/{interval}_vep_out.txt",
#    log:
#        os.path.join(
#            config["workdir"],
#            "logs/run_vep_invnorm/hdrr/dist_angle/{interval}/15/5000/0.8/invnorm.log"
#        ),
#    params:
#        species = lambda wildcards: config["ref"]["species"],
#        assembly = lambda wildcards: config["ref"]["build"],
#    resources:
#        mem_mb = 1000
#    container:
#        config["ensembl_vep_104"]
#    shell:
#        """
#        vep \
#            --input_file {input} \
#            --output_file {output.out} \
#            --species {params.species} \
#            --assembly {params.assembly} \
#            --database \
#            --force_overwrite \
#                2> {log}
#        """

rule sig_snps_boxplots:
    input:
        bed = "/hps/nobackup/birney/users/ian/MIKK_HMM/beds/F2/hdrr/5000/0.8.bed",
        bim = "/hps/nobackup/birney/users/ian/MIKK_HMM/beds/F2/hdrr/5000/0.8.bim",
        fam = "/hps/nobackup/birney/users/ian/MIKK_HMM/beds/F2/hdrr/5000/0.8.fam",
        line_cols = rules.get_line_ranks_and_colours.output.csv,
        samples_file = config["F2_samples_file"],
        sigs = os.path.join(
            config["workdir"],
            "sig_snps/hdrr/dist_angle/{interval}/15/5000/0.8/invnorm/sigs.csv"
        ),
        phenos = expand(os.path.join(
            config["workdir"],
            "phens_sf/true/dist_angle/{{interval}}/15/{dge_sge}/invnorm/{assay}/{state}.phen",
                ),
                dge_sge = ["dge", "sge"],
                assay = ["open_field", "novel_object"],
                state = list(range(1,16))
        ),
    output:
        greedy = os.path.join(
            config["workdir"],
            "sig_snps/hdrr/dist_angle/{interval}/15/5000/0.8/invnorm/sigs_greedy_filt.csv"
        ),
    log:
        os.path.join(
            config["workdir"],
            "logs/sig_snps_boxplots/dist_angle/{interval}/15/5000/0.8/invnorm.log"
        ),
    params:
        out_dir = lambda wildcards: os.path.join(
            "book/figs/sig_snps_boxplots/dist_angle",
            wildcards.interval,
            "15/invnorm"
        ),
        n_states = 15
    resources:
        mem_mb = 60000
    container:
        config["R_4.2.0"]
    script:
        "../scripts/sig_snps_boxplots.R"

