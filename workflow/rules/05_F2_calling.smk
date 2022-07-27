# Get read counts supporting each F0-homozygous-divergent allele
rule bam_readcount_F2:
    input:
        bam = rules.merge_bams.output.bam,
        index = rules.samtools_index.output,
        sites_file = rules.extract_homo_div_snps.output.sites,
        ref = os.path.join(
            config["workdir"],
            "refs/hdrr.fasta"
        ),
    output:
        os.path.join(
            config["workdir"],
            "dp4s/F2/hdrr/{sample}_{pat}_{mat}.dp4.txt"
        ),
    log:
        os.path.join(
            config["workdir"],
            "logs/bam_readcount_F2/hdrr/{sample}_{pat}_{mat}.log"
        ),
    resources:
        mem_mb = 200
    container:
        config["bam-readcount"]
    shell:
        """
        bam-readcount \
            -l {input.sites_file} \
            -f {input.ref} \
            {input.bam} | \
            cut -f 1,15,28,41,54,67 -d ":" | sed 's/=//g' | sed 's/\\t:/\\t/g' | sed 's/:/\\t/g' \
                > {output} 2> {log}
        """

# Make dpAB files
rule make_dp_AB_F2:
    input:
        dp4 = rules.bam_readcount_F2.output,
        sites_file = rules.extract_homo_div_snps.output.sites,
    output:
        os.path.join(
            config["workdir"],
            "dpABs/F2/hdrr/{sample}_{pat}_{mat}.txt"
        ),
    log:
        os.path.join(
            config["workdir"],
            "logs/make_dp_AB_F2/hdrr/{sample}_{pat}_{mat}.log"
        ),
    resources:
        mem_mb = 2000
    script:
        "../scripts/make_dp_AB.py"

# Create HMM inputs
rule make_hmm_input:
    input:
        expand(rules.make_dp_AB_F2.output,
            zip,
            sample = SAMPLES_F2_ZIP,
            pat = PAT_ZIP,
            mat = MAT_ZIP
        ),
    output:
        csv = os.path.join(
            config["workdir"],
            "hmm_in/F2/hdrr/{bin_length}.csv"
        ),
    log:
        os.path.join(
            config["workdir"],
            "logs/make_hmm_input/hdrr/{bin_length}.log"
        ),
    params:
        bin_length = "{bin_length}",
    resources:
        mem_mb = 7000
    container:
        config["R_4.2.0"]
    script:
        "../scripts/make_hmm_input.R"

# Run HMM to genotype F2
rule true_hmmlearn:
    input:
        rules.make_hmm_input.output,
    output:
        csv = os.path.join(
            config["workdir"],
            "hmm_out/F2/hdrr/hmmlearn_true/{bin_length}/{cov}.csv"
        ),
        pck = os.path.join(
            config["workdir"],
            "hmm_out/F2/hdrr/hmmlearn_true/{bin_length}/{cov}.pickle"
        ),
    log:
        os.path.join(
            config["workdir"],
            "logs/true_hmmlearn/hdrr/{bin_length}/{cov}.log"
        ),
    params:
        cov = "{cov}",
        #low_cov_samples = lambda wildcards: config["low_cov_samples"]
    resources:
        mem_mb = 20000
    container:
        config["hmmlearn"]
    script:
        "../scripts/true_hmmlearn.py"

# Split genotyped F2 by sample
rule split_HMM_genotyped_F2:
    input:
        hmm = rules.true_hmmlearn.output.csv,
    output:
        expand(os.path.join(
            config["workdir"],
            "hmm_out/F2/hdrr/hmmlearn_split/{{bin_length}}/{{cov}}/{sample}_{pat}_{mat}.csv"
            ),
                zip,
                sample = SAMPLES_F2_ZIP,
                pat = PAT_ZIP,
                mat = MAT_ZIP
        )
    log:
        os.path.join(
            config["workdir"],
            "logs/split_HMM_genotyped_F2/hdrr/{bin_length}/{cov}.log"
        ),
    resources:
        mem_mb = 50000
    container:
        config["R_4.2.0"]
    script:
        "../scripts/split_HMM_genotyped_F2.R"

# Take HMM genotypes for F2 and map them back to the genotypes of the F0
rule impute_F2_genos:
    input:
        F2 = os.path.join(
            config["workdir"],
            "hmm_out/F2/hdrr/hmmlearn_split/{bin_length}/{cov}/{sample}_{pat}_{mat}.csv"
            ),
        F0 = rules.extract_homo_div_snps.output.sites,
    output:
        nt = os.path.join(
            config["workdir"],
            "F2_with_genos/hdrr/{bin_length}/{cov}/{sample}_{pat}_{mat}.csv"
        ),
        ab = os.path.join(
            config["workdir"],
            "F2_with_genos_AB/hdrr/{bin_length}/{cov}/{sample}_{pat}_{mat}.csv"
        ),
    log:
        os.path.join(
            config["workdir"],
            "logs/impute_F2_genos/hdrr/{bin_length}/{cov}/{sample}_{pat}_{mat}.log"
        ),
    params:
        bin_length = "{bin_length}"
    resources:
        mem_mb = 3000
    container:
        config["R_4.2.0"]
    script:
        "../scripts/impute_F2_genos.R"

