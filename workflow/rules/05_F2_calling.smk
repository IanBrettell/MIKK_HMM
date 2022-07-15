# Get read counts supporting each F0-homozygous-divergent allele
rule bam_readcount_F2:
    input:
        bam = rules.merge_bams.output.bam,
        index = rules.samtools_index.output,
        sites_file = rules.extract_homo_div_snps.output.sites,
        ref = rules.get_genome.output,
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
