#rule get_genome:
#    output:
#        os.path.join(
#            config["workdir"],
#            "refs/hdrr.fasta"
#        ),
#    log:
#        os.path.join(
#            config["workdir"],
#            "logs/get_genome/hdrr.log"
#        ),
#    params:
#        species=lambda wildcards: config["ref"]["species"],
#        datatype="dna",
#        build=lambda wildcards: config["ref"]["build"],
#        release=lambda wildcards: config["ref"]["release"],
#    resources:
#        mem_mb = 1000
#    wrapper:
#        "v1.7.0/bio/reference/ensembl-sequence"

rule genome_faidx:
    input:
        os.path.join(
            config["workdir"],
            "refs/hdrr.fasta"
        ),
    output:
        os.path.join(
            config["workdir"],
            "refs/hdrr.fasta.fai"
        )
    log:
        os.path.join(
            config["workdir"],
            "logs/genome_faidx/hdrr.log"
        ),
    container:
        config["samtools"]
    resources:
        mem_mb = 1000
    shell:
        """
        samtools faidx {input[0]}  > {output[0]} \
            2> {log}
        """

rule genome_dict:
    input:
        os.path.join(
            config["workdir"],
            "refs/hdrr.fasta"
        ),
    output:
        os.path.join(
            config["workdir"],
            "refs/hdrr.dict"
        ),
    log:
        os.path.join(
            config["workdir"],
            "logs/genome_dict/hdrr.log"
        ),
    container:
        config["samtools"]
    resources:
        mem_mb = 5000
    shell:
        """
        samtools dict {input} > {output} \
            2> {log}
        """

rule bwa_mem2_index:
    input:
        os.path.join(
            config["workdir"],
            "refs/hdrr.fasta"
        ),
    output:
        multiext(
            os.path.join(
                config["workdir"],
                "refs/hdrr.fasta"),
            ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
    log:
        os.path.join(
            config["workdir"],
            "logs/bwa_mem2_index/hdrr.log"
        ),    
    resources:
        mem_mb=50000,
    container:
        config["bwa-mem2"]
    shell:
        """
        bwa-mem2 index {input[0]} \
            2> {log}
        """

rule get_chrom_lengths:
    input:
        os.path.join(
            config["workdir"],
            "refs/hdrr.fasta"
        ),
    output:
        csv = "config/hdrr_chrom_lengths.csv"
    log:
        os.path.join(
            config["workdir"],
            "logs/get_chrom_lengths/hdrr.log"
        ),
    container:
        config["bash"]
    resources:
        mem_mb = 200
    shell:
        """
        grep "^>" {input} | cut -f1 -d' ' | sed 's/>//g' > tmp1.txt ;
        grep "^>" {input} | cut -f3 -d' ' | cut -f5 -d':' > tmp2.txt ;
        paste -d',' tmp1.txt tmp2.txt > {output.csv} ;
        rm tmp1.txt tmp2.txt \
            2> {log}
        """
    