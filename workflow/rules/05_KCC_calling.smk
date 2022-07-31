# Call haplotypes
rule haplotype_caller:
    input:
        bam=rules.merge_bams.output,
        bai=rules.samtools_index.output,
        ref= os.path.join(
            config["workdir"],
            "refs/hdrr.fasta"
        ),
        ref_index = rules.genome_faidx.output,
        ref_dict = rules.genome_dict.output
    output:
        gvcf=os.path.join(
            config["workdir"],
            "vcfs/KCC/hdrr/gvcfs/{sample}/{contig}.g.vcf"
        ),
    log:
        os.path.join(
            config["workdir"],
            "logs/haplotype_caller/hdrr/{sample}_{contig}.log"
        )
    params:
        extra = lambda wildcards: "-L " + wildcards.contig + " --tmp-dir " + config["tmpdir"]
    resources:
        java_mem_gb=4,
        mem_mb=10000
    container:
        config["gatk"]
    shell:
        """
        gatk --java-options \"-Xmx{resources.java_mem_gb}G\" \
            HaplotypeCaller \
                {params.extra} \
                -R {input.ref} \
                -I {input.bam} \
                -ERC GVCF \
                -O {output.gvcf} \
                    > {log} 2>&1
        """

# Combine all samples into 
rule combine_calls:
    input:
        ref= os.path.join(
            config["workdir"],
            "refs/hdrr.fasta"
        ),
        gvcfs=expand(os.path.join(
                config["workdir"],
                "vcfs/KCC/hdrr/gvcfs/{sample}/{{contig}}.g.vcf"),
                    sample=SAMPLES_KCC
        ),
    output:
        os.path.join(
            config["workdir"],
            "vcfs/KCC/hdrr/combined/all.{contig}.g.vcf.gz"
        ),
    log:
        os.path.join(
            config["workdir"], 
            "logs/combine_calls/hdrr/{contig}.log"
        ),
    params:
        in_files = lambda wildcards, input: " -V ".join(input.gvcfs)
    resources:
        java_mem_gb=1,
        mem_mb = 5000
    container:
        config["gatk"]
    shell:
        """
        gatk --java-options \"-Xmx{resources.java_mem_gb}G\"\
            CombineGVCFs \
                -V {params.in_files} \
                -R {input.ref} \
                -O {output[0]} \
                    > {log} 2>&1
        """

rule genotype_variants:
    input:
        ref= os.path.join(
            config["workdir"],
            "refs/hdrr.fasta"
        ),
        gvcf=rules.combine_calls.output,
    output:
        os.path.join(
            config["workdir"], 
            "vcfs/KCC/hdrr/genotyped/all.{contig}.vcf.gz"
        ),
    log:
        os.path.join(
            config["workdir"], 
            "logs/genotype_variants/hdrr/{contig}.log"
        ),
    resources:
        java_mem_gb=1,
        mem_mb = 5000
    container:
        config["gatk"]
    shell:
        """
        gatk --java-options \"-Xmx{resources.java_mem_gb}G\" \
            GenotypeGVCFs \
                -V {input.gvcf} \
                -R {input.ref} \
                -O {output[0]} \
                    > {log} 2>&1
        """

rule merge_variants:
    input:
        vcfs=expand(os.path.join(
                config["workdir"], 
                "vcfs/KCC/hdrr/genotyped/all.{contig}.vcf.gz"),
                    contig=get_contigs()
        ),
    output:
        os.path.join(
            config["workdir"], 
            "vcfs/KCC/hdrr/final/all.vcf.gz"
        ),
    log:
        os.path.join(
            config["workdir"], 
            "logs/merge_variants/hdrr.log"
        ),
    params:
        in_files = lambda wildcards, input: " ".join("INPUT={}".format(f) for f in input.vcfs)
    resources:
        mem_mb = 5000
    container:
        config["picard"]
    shell:
        """
        picard MergeVcfs \
            {params.in_files} \
            OUTPUT={output[0]} \
                &> {log}
        """
