# Generate genetic relationship matrix
rule make_grm:
    input:
        bed = rules.create_bed.output.bed,
    output:
        os.path.join(
            config["workdir"],
            "grms/hdrr/{bin_length}/{cov}.grm.bin"
        ),     
    log:
        os.path.join(
            config["workdir"],
            "logs/make_grm/hdrr/{bin_length}/{cov}.log"
        ),
    params:
        in_pref = lambda wildcards, input: input.bed.replace(".bed", ""),
        out_pref = lambda wildcards, output: output[0].replace(".grm.bin", ""),    
    resources:
        mem_mb = 5000,
        threads = 1
    container:
        config["GCTA"]
    shell:
        """
        gcta64 \
            --make-grm \
            --bfile {params.in_pref} \
            --out {params.out_pref} \
            --autosome-num 24 \
            --thread-num {resources.threads} \
                2> {log}
        """

# Generate genetic relationship matrix for each chromosome
rule make_grm_per_chrom:
    input:
        bed = rules.create_bed.output.bed,
    output:
        os.path.join(
            config["workdir"],
            "grms_per_chr/hdrr/{bin_length}/{cov}/{contig}.grm.bin"
        ),     
    log:
        os.path.join(
            config["workdir"],
            "logs/make_grm_per_chrom/hdrr/{bin_length}/{cov}/{contig}.log"
        ),
    params:
        in_pref = lambda wildcards, input: input.bed.replace(".bed", ""),
        out_pref = lambda wildcards, output: output[0].replace(".grm.bin", ""),
        contig = "{contig}"   
    resources:
        mem_mb = 5000,
        threads = 1
    container:
        config["GCTA"]
    shell:
        """
        gcta64 \
            --make-grm \
            --chr {params.contig} \
            --bfile {params.in_pref} \
            --out {params.out_pref} \
            --autosome-num 24 \
            --thread-num {resources.threads} \
                2> {log}
        """

rule make_grm_inbred:
    input:
        bed = rules.create_bed.output.bed,
    output:
        os.path.join(
            config["workdir"],
            "grms_inbred/hdrr/{bin_length}/{cov}.grm.bin"
        ),     
    log:
        os.path.join(
            config["workdir"],
            "logs/make_grm_inbred/hdrr/{bin_length}/{cov}.log"
        ),
    params:
        in_pref = lambda wildcards, input: input.bed.replace(".bed", ""),
        out_pref = lambda wildcards, output: output[0].replace(".grm.bin", ""),    
    resources:
        mem_mb = 80000,
        threads = 1
    container:
        config["GCTA"]
    shell:
        """
        gcta64 \
            --make-grm-inbred \
            --bfile {params.in_pref} \
            --out {params.out_pref} \
            --autosome-num 24 \
            --thread-num {resources.threads} \
                2> {log}
        """

rule make_grm_inbred_per_chrom:
    input:
        bed = rules.create_bed.output.bed,
    output:
        os.path.join(
            config["workdir"],
            "grms_inbred_per_chr/hdrr/{bin_length}/{cov}/{contig}.grm.bin"
        ),     
    log:
        os.path.join(
            config["workdir"],
            "logs/make_grm_inbred_per_chrom/hdrr/{bin_length}/{cov}/{contig}.log"
        ),
    params:
        in_pref = lambda wildcards, input: input.bed.replace(".bed", ""),
        out_pref = lambda wildcards, output: output[0].replace(".grm.bin", ""),
        contig = "{contig}" 
    resources:
        mem_mb = 80000,
        threads = 1
    container:
        config["GCTA"]
    shell:
        """
        gcta64 \
            --make-grm-inbred \
            --chr {params.contig} \
            --bfile {params.in_pref} \
            --out {params.out_pref} \
            --autosome-num 24 \
            --thread-num {resources.threads} \
                2> {log}
        """

rule plink_grm:
    input:
        rules.create_bed.output.bed
    output:
        dist = os.path.join(
            config["workdir"],
            "grm_plink/F2/hdrr/{bin_length}/{cov}/F2/all.grm.bin"
        ),
    log:
        os.path.join(
            config["workdir"],
            "logs/plink_grm/hdrr/{bin_length}/{cov}/all.log"
        ),
    params:
        in_pref = lambda wildcards, input: input[0].replace(".bed", ""),
        out_pref = lambda wildcards, output: output.dist.replace(".dist", ""),
        contig = "{contig}"
    resources:
        mem_mb = 5000
    container:
        config["plink1.9"]
    shell:
        """
        plink1.9 \
            --bfile {params.in_pref} \
            --make-grm-bin \
            --no-fid \
            --no-parents \
            --no-sex \
            --no-pheno \
            --chr-set 24 no-xy no-mt \
            --out {params.out_pref} \
                2> {log}
        """

rule plink_grm_loco:
    input:
        rules.create_bed.output.bed
    output:
        dist = os.path.join(
            config["workdir"],
            "grms_loco/F2/hdrr/{bin_length}/{cov}/F2/{contig}.grm.bin"
        ),
    log:
        os.path.join(
            config["workdir"],
            "logs/plink_grm_loco/hdrr/{bin_length}/{cov}/{contig}.log"
        ),
    params:
        in_pref = lambda wildcards, input: input[0].replace(".bed", ""),
        out_pref = lambda wildcards, output: output.dist.replace(".dist", ""),
        contig = "{contig}"
    resources:
        mem_mb = 5000
    container:
        config["plink1.9"]
    shell:
        """
        plink1.9 \
            --bfile {params.in_pref} \
            --make-grm-bin \
            --not-chr {params.contig} \
            --no-fid \
            --no-parents \
            --no-sex \
            --no-pheno \
            --chr-set 24 no-xy no-mt \
            --out {params.out_pref} \
                2> {log}
        """

rule plink_distance:
    input:
        rules.create_bed.output.bed
    output:
        dist = os.path.join(
            config["workdir"],
            "relationships/F2/hdrr/{bin_length}/{cov}/F2.dist"
        ),
    log:
        os.path.join(
            config["workdir"],
            "logs/plink_distance/hdrr/{bin_length}/{cov}.log"
        ),
    params:
        in_pref = lambda wildcards, input: input[0].replace(".bed", ""),
        out_pref = lambda wildcards, output: output.dist.replace(".dist", ""),
    resources:
        mem_mb = 5000
    container:
        config["plink1.9"]
    shell:
        """
        plink1.9 \
            --bfile {params.in_pref} \
            --distance square \
            --no-fid \
            --no-parents \
            --no-sex \
            --no-pheno \
            --chr-set 24 no-xy no-mt \
            --out {params.out_pref} \
                2> {log}
        """

rule gcta_pca:
    input:
        rules.make_grm_inbred.output,
    output:
        os.path.join(
            config["workdir"],
            "gcta_pca/hdrr/{bin_length}/{cov}.eigenval"
        ),     
    log:
        os.path.join(
            config["workdir"],
            "logs/gcta_pca/hdrr/{bin_length}/{cov}.log"
        ),
    params:
        in_pref = lambda wildcards, input: input[0].replace(".grm.bin", ""),
        out_pref = lambda wildcards, output: output[0].replace(".eigenval", ""),    
    resources:
        mem_mb = 10000,
        threads = 1
    container:
        config["GCTA"]
    shell:
        """
        gcta64 \
            --pca \
            --grm {params.in_pref} \
            --out {params.out_pref} \
            --autosome-num 24 \
            --thread-num {resources.threads} \
                2> {log}
        """

#####################
# No missing
#####################

# Generate genetic relationship matrix
rule make_grm_no_miss:
    input:
        bed = rules.create_bed_no_miss.output.bed,
    output:
        os.path.join(
            config["workdir"],
            "grms_no_miss/hdrr/{bin_length}/{cov}.grm.bin"
        ),     
    log:
        os.path.join(
            config["workdir"],
            "logs/make_grm_no_miss/hdrr/{bin_length}/{cov}.log"
        ),
    params:
        in_pref = lambda wildcards, input: input.bed.replace(".bed", ""),
        out_pref = lambda wildcards, output: output[0].replace(".grm.bin", ""),    
    resources:
        mem_mb = 5000,
        threads = 1
    container:
        config["GCTA"]
    shell:
        """
        gcta64 \
            --make-grm \
            --bfile {params.in_pref} \
            --out {params.out_pref} \
            --autosome-num 24 \
            --thread-num {resources.threads} \
                2> {log}
        """

rule make_grm_inbred_no_miss:
    input:
        bed = rules.create_bed_no_miss.output.bed,
    output:
        os.path.join(
            config["workdir"],
            "grms_inbred_no_miss/hdrr/{bin_length}/{cov}.grm.bin"
        ),     
    log:
        os.path.join(
            config["workdir"],
            "logs/make_grm_inbred_no_miss/hdrr/{bin_length}/{cov}.log"
        ),
    params:
        in_pref = lambda wildcards, input: input.bed.replace(".bed", ""),
        out_pref = lambda wildcards, output: output[0].replace(".grm.bin", ""),    
    resources:
        mem_mb = 80000,
        threads = 1
    container:
        config["GCTA"]
    shell:
        """
        gcta64 \
            --make-grm-inbred \
            --bfile {params.in_pref} \
            --out {params.out_pref} \
            --autosome-num 24 \
            --thread-num {resources.threads} \
                2> {log}
        """
