# Create .ped and .map files with format
rule create_ped:
    input:
        expand(os.path.join(
            config["workdir"],
            "F2_with_genos/hdrr/{{bin_length}}/{{cov}}/{sample}_{pat}_{mat}.csv"
            ),
                zip,
                sample = SAMPLES_F2_ZIP,
                pat = PAT_ZIP,
                mat = MAT_ZIP
        ),
    output:
        ped = os.path.join(
            config["workdir"],
            "peds/F2/hdrr/{bin_length}/{cov}.ped"
        ),
        map = os.path.join(
            config["workdir"],
            "peds/F2/hdrr/{bin_length}/{cov}.map"
        ),
    log:
        os.path.join(
            config["workdir"],
            "logs/create_ped/{bin_length}/{cov}.log"
        ),
    resources:
        mem_mb = 80000,
    container:
        # requires tidyr >= v1.2
        config["tidyverse_4.1.3"]
    script:
        "../scripts/create_ped.R"

# Convert .ped to .bed
rule create_bed:
    input:
        rules.create_ped.output.ped
    output:
        bed = os.path.join(
            config["workdir"],
            "beds/F2/hdrr/{bin_length}/{cov}.bed"
        ),
    log:
        os.path.join(
            config["workdir"],
            "logs/create_bed/hdrr/{bin_length}/{cov}.log"
        ),
    params:
        in_pref = lambda wildcards, input: input[0].replace(".ped", ""),
        out_pref = lambda wildcards, output: output.bed.replace(".bed", ""),
    resources:
        mem_mb = 1000
    container:
        config["plink1.9"]
    shell:
        """
        plink1.9 \
            --make-bed \
            --no-fid \
            --no-parents \
            --no-sex \
            --no-pheno \
            --chr-set 24 no-xy no-mt \
            --file {params.in_pref} \
            --out {params.out_pref} \
                2> {log}
        """

# Create .phen files
# as specified here: https://gcta.freeforums.net/thread/247/greml-estimating-variance-explained-snps
rule create_phen:
    input:
        ped = expand(rules.create_ped.output,
            bin_length = 5000,
            cov = 0.8
        ),
        phenos = os.path.join(
            config["workdir"],
            "state_freq_F2/dist_angle/0.05/15/{dge_sge}/{transformation}/{state}.csv"
            ),
        samples_file = config["F2_samples_file"]
    output:
        os.path.join(
            config["workdir"],
            "phens/true/{dge_sge}/{transformation}/{assay}/{state}.phen"
        ),
    log:
        os.path.join(
            config["workdir"],
            "logs/create_phen/{dge_sge}/{transformation}/{assay}/{state}.log"
        ),
    params:
        assay = "{assay}"        
    resources:
        mem_mb = 20000,
    container:
        # requires tidyr >= v1.2
        config["tidyverse_4.1.3"]
    script:
        "../scripts/create_phen.R"

rule permute_phen:
    input:
        rules.create_phen.output,
    output:
        os.path.join(
            config["workdir"],
            "phens/permuted/{dge_sge}/{transformation}/{assay}/{state}/{seed}.phen"
        ),
    log:
        os.path.join(
            config["workdir"],
            "logs/permute_phen/{dge_sge}/{transformation}/{assay}/{state}/{seed}.log"
        ),
    params:
        phenotype = "state_freq",
        seed = "{seed}"
    resources:
        mem_mb = 500,
    container:
        # requires tidyr >= v1.2
        config["tidyverse_4.1.3"]
    script:
        "../scripts/permute_phen.R"
