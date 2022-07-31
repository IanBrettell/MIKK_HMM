# Get samples with no calls for a given chromosome
rule get_samples_no_calls:
    input:
        expand(os.path.join(
            config["workdir"],
            "F2_with_genos_AB_chr/hdrr/{{bin_length}}/{{cov}}/{sample}_{pat}_{mat}/{{contig}}.csv"
            ),
                zip,
                sample = SAMPLES_F2_ZIP,
                pat = PAT_ZIP,
                mat = MAT_ZIP
        ),
    output:
        csv = "config/samples_no_calls/{bin_length}/{cov}/{contig}.csv"
    log:
        os.path.join(
            config["workdir"],
            "logs/get_samples_no_calls/{bin_length}/{cov}/{contig}.log"
        ),
    params:
        contig = "{contig}"
    resources:
        mem_mb = 30000,
    container:
        config["R_4.2.0"]
    script:
        "../scripts/get_samples_no_calls.R"        

# Get full list of samples with no calls on any chromosome
rule get_samples_no_calls_all_chr:
    input:
        expand("config/samples_no_calls/{{bin_length}}/{{cov}}/{contig}.csv",
                contig = list(range(1,25))
        ),
    output:
        csv = "config/samples_no_calls_all/{bin_length}/{cov}.csv"
    log:
        os.path.join(
            config["workdir"],
            "logs/get_samples_no_calls/{bin_length}/{cov}.log"
        ),
    resources:
        mem_mb = 30000,
    container:
        # requires tidyr >= v1.2
        config["R_4.2.0"]
    script:
        "../scripts/get_samples_no_calls_all_chr.R"        

# Get list of samples do exclude due to no calls in a given chromosome
S_NO_CALLS = pd.read_csv("config/samples_no_calls_all/5000/0.8.csv")
S_NO_CALLS = S_NO_CALLS['SAMPLE'].tolist()

# Create .ped and .map files with format
rule create_ped_contigs:
    input:
        expand(os.path.join(
            config["workdir"],
            "F2_with_genos_AB_chr/hdrr/{{bin_length}}/{{cov}}/{sample}_{pat}_{mat}/{{contig}}.csv"
            ),
                zip,
                sample = SAMPLES_F2_ZIP,
                pat = PAT_ZIP,
                mat = MAT_ZIP
        ),
    output:
        ped = os.path.join(
            config["workdir"],
            "peds_contigs/F2/hdrr/{bin_length}/{cov}/{contig}.ped"
        ),
        map = os.path.join(
            config["workdir"],
            "peds_contigs/F2/hdrr/{bin_length}/{cov}/{contig}.map"
        ),
    log:
        os.path.join(
            config["workdir"],
            "logs/create_ped/{bin_length}/{cov}/{contig}.log"
        ),
    resources:
        mem_mb = 30000,
    container:
        # requires tidyr >= v1.2
        config["tidyverse_4.1.3"]
    script:
        "../scripts/create_ped.R"

# Remove sample columns
rule remove_sample_col:
    input:
        rules.create_ped_contigs.output.ped
    output:
        ped = os.path.join(
            config["workdir"],
            "peds_contigs_no_sample/F2/hdrr/{bin_length}/{cov}/{contig}.ped"
        ),
        ids = os.path.join(
            config["workdir"],
            "peds_contigs_sample_ids/F2/hdrr/{bin_length}/{cov}/{contig}.ids"
        ),
    log:
        os.path.join(
            config["workdir"],
            "logs/remove_sample_col/{bin_length}/{cov}/{contig}.log"
        ),
    params:
        contig = "{contig}"
    resources:
        mem_mb = 5000,
    shell:
        """
        cut -f1 {input} > {output.ids} && \
        cut -f2- {input} > {output.ped} \
            2> {log}
        """

# Combine contig .peds into single .ped
rule combine_peds:
    input:
        ids = os.path.join(
            config["workdir"],
            "peds_contigs_sample_ids/F2/hdrr/{bin_length}/{cov}/1.ids"
        ),
        peds = expand(os.path.join(
            config["workdir"],
            "peds_contigs_no_sample/F2/hdrr/{{bin_length}}/{{cov}}/{contig}.ped"
            ),
                contig = list(range(1, 25))
        ),
        maps = expand(os.path.join(
            config["workdir"],
            "peds_contigs/F2/hdrr/{{bin_length}}/{{cov}}/{contig}.map"
            ),
                contig = list(range(1, 25))
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
            "logs/combine_peds/{bin_length}/{cov}.log"
        ),
    resources:
        mem_mb = 50000,
    shell:
        """
        paste -d$'\\t' {input.ids} {input.peds} > {output.ped} && \
        cat {input.maps} > {output.map} \
            2> {log}
        """

# Convert .ped to .bed
rule create_bed:
    input:
        rules.combine_peds.output.ped
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
        mem_mb = 5000
    container:
        config["plink1.9"]
    shell:
        """
        plink1.9 \
            --make-bed \
            --geno 0 \
            --no-fid \
            --no-parents \
            --no-sex \
            --no-pheno \
            --chr-set 24 no-xy no-mt \
            --file {params.in_pref} \
            --out {params.out_pref} \
                2> {log}
        """

rule test_grm:
    input:
        expand(rules.impute_F2_genos.output.ab,
            zip,
            bin_length = 5000,
            cov = 0.8,
            contig = 10,
            sample = [2, 1, 3, 18, 17, 20],
            pat = ["38-2", "38-2", "38-2", "8-2", "8-2", "8-2"],
            mat = ["21-2", "21-2", "21-2", "40-1", "40-1", "40-1"]
        )
    output:
        toy_ped = os.path.join(
            config["workdir"],
            "grm_test/F2/hdrr/{bin_length}/{cov}.bed"
        ),
    log:
        os.path.join(
            config["workdir"],
            "logs/create_bed/hdrr/{bin_length}/{cov}.log"
        ),
    resources:
        mem_mb = 30000,
    script:
        "../scripts/test_grm.R" 

#####################
# No missing
#####################

# Create .ped and .map files with format
## This time, remove all SNPs with a single missing genotype
rule create_ped_contigs_no_miss:
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
            "peds_contigs_no_miss/F2/hdrr/{bin_length}/{cov}/{contig}.ped"
        ),
        map = os.path.join(
            config["workdir"],
            "peds_contigs_no_miss/F2/hdrr/{bin_length}/{cov}/{contig}.map"
        ),
    log:
        os.path.join(
            config["workdir"],
            "logs/create_ped_contigs_no_miss/{bin_length}/{cov}/{contig}.log"
        ),
    params:
        contig = "{contig}"
    resources:
        mem_mb = 30000,
    container:
        # requires tidyr >= v1.2
        config["tidyverse_4.1.3"]
    script:
        "../scripts/create_ped_no_miss.R"

# Remove sample columns
rule remove_sample_col_no_miss:
    input:
        rules.create_ped_contigs_no_miss.output.ped
    output:
        ped = os.path.join(
            config["workdir"],
            "peds_contigs_no_sample_no_miss/F2/hdrr/{bin_length}/{cov}/{contig}.ped"
        ),
        ids = os.path.join(
            config["workdir"],
            "peds_contigs_sample_ids_no_miss/F2/hdrr/{bin_length}/{cov}/{contig}.ids"
        ),
    log:
        os.path.join(
            config["workdir"],
            "logs/remove_sample_col_no_miss/{bin_length}/{cov}/{contig}.log"
        ),
    params:
        contig = "{contig}"
    resources:
        mem_mb = 5000,
    shell:
        """
        cut -f1 {input} > {output.ids} && \
        cut -f2- {input} > {output.ped} \
            2> {log}
        """

# Combine contig .peds into single .ped
rule combine_peds_no_miss:
    input:
        ids = os.path.join(
            config["workdir"],
            "peds_contigs_sample_ids_no_miss/F2/hdrr/{bin_length}/{cov}/1.ids"
        ),
        peds = expand(os.path.join(
            config["workdir"],
            "peds_contigs_no_sample_no_miss/F2/hdrr/{{bin_length}}/{{cov}}/{contig}.ped"
            ),
                contig = list(range(1, 25))
        ),
        maps = expand(os.path.join(
            config["workdir"],
            "peds_contigs_no_miss/F2/hdrr/{{bin_length}}/{{cov}}/{contig}.map"
            ),
                contig = list(range(1, 25))
        ),
    output: 
        ped = os.path.join(
            config["workdir"],
            "peds_no_miss/F2/hdrr/{bin_length}/{cov}.ped"
        ),
        map = os.path.join(
            config["workdir"],
            "peds_no_miss/F2/hdrr/{bin_length}/{cov}.map"
        ),
    log:
        os.path.join(
            config["workdir"],
            "logs/combine_peds_no_miss/{bin_length}/{cov}.log"
        ),
    resources:
        mem_mb = 50000,
    shell:
        """
        paste -d$'\\t' {input.ids} {input.peds} > {output.ped} && \
        cat {input.maps} > {output.map} \
            2> {log}
        """

# Convert .ped to .bed
rule create_bed_no_miss:
    input:
        rules.combine_peds_no_miss.output.ped
    output:
        bed = os.path.join(
            config["workdir"],
            "beds_no_miss/F2/hdrr/{bin_length}/{cov}.bed"
        ),
    log:
        os.path.join(
            config["workdir"],
            "logs/create_bed_no_miss/hdrr/{bin_length}/{cov}.log"
        ),
    params:
        in_pref = lambda wildcards, input: input[0].replace(".ped", ""),
        out_pref = lambda wildcards, output: output.bed.replace(".bed", ""),
    resources:
        mem_mb = 5000
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

