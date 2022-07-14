# Create samples file with sequence paths
# NOTE: this needs to be run alone, before any of the following rules
rule create_seq_sample_file:
    input:
        in_dir = os.path.join(
            config["sequence_dir"],
            "{run}",
            "aspera_{run}"
        )
    output:
        out = "config/seq_sample_file_{run}.txt"
    log:
        os.path.join(
            config["workdir"], 
            "logs/create_seq_sample_file/{run}.log"
        )
    container:
        config["R_4.2.0"]
    resources:
        mem_mb = 500
    script:
        "../scripts/create_seq_sample_file.R"

# Get list of samples
# NOTE: this is dependent on both (all?) runs having the same samples
SAMPLES = pd.read_csv("config/seq_sample_file_AAAMKYVHV.txt",
                      sep = '\t')['sample'].values.tolist()

# Get paired sequence files for each run/sample combination
def get_seq_paths(wildcards):
    df = pd.read_csv("config/seq_sample_file_" + wildcards.run + ".txt",
                     sep = '\t')
    target_row = df.loc[df['sample'] == wildcards.sample]
    fq1 = target_row['fq1'].values[0]
    fq2 = target_row['fq2'].values[0]
    return [fq1, fq2]

# Align sequences
rule bwa_mem2_mem:
    input:
        reads=get_seq_paths,
        idx=rules.bwa_mem2_index.output,
        ref = rules.get_genome.output,
    output:
        os.path.join(
            config["workdir"],
            "sams/hdrr/bwamem2/mapped/{run}/{sample}.sam"
        ),
    log:
        os.path.join(
            config["workdir"], 
            "logs/bwa_mem2_mem/hdrr/{run}/{sample}.log"
        ),
    params:
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
    container:
        config["bwa-mem2"]
    resources:
        mem_mb = 10000
    threads:
        8
    shell:
        """
        bwa-mem2 mem \
            -t {threads} \
            {params.extra} \
            {input.ref} \
            {input.reads} \
                > {output} \
                    2> {log}
        """

rule sort_sam:
    input:
        rules.bwa_mem2_mem.output,
    output:
        os.path.join(
            config["workdir"],
            "bams/bwamem2/sorted/{run}/{sample}.bam"
        ),
    log:
        os.path.join(
            config["workdir"], 
            "logs/sort_sam/{run}/{sample}.log"
        ),
    params:
        sort_order="coordinate",
        extra=lambda wildcards: "VALIDATION_STRINGENCY=LENIENT tmpdir=" + config["tmpdir"],
    container:
        config["picard"]
    resources:
        java_mem_mb = 1024,
        mem_mb = 2000
    shell:
        """
        picard SortSam \
            -Xmx{resources.java_mem_mb}M \
            {params.extra} \
            INPUT={input[0]} \
            OUTPUT={output[0]} \
            SORT_ORDER={params.sort_order} \
                2> {log}
        """
rule mark_duplicates:
    input:
        rules.sort_sam.output,
    output:
        bam=os.path.join(
            config["workdir"], 
            "bams/hdrr/bwamem2/marked/{run}/{sample}.bam"
        ),
        metrics=os.path.join(
            config["workdir"], 
            "bams/hdrr/bwamem2/marked/{run}/{sample}.metrics.txt"
        ),
    log:
        os.path.join(
            config["workdir"],
            "logs/mark_duplicates/hdrr/{run}/{sample}.log"
        ),
    params:
        extra = lambda wildcards: "REMOVE_DUPLICATES=true tmpdir=" + config["tmpdir"]
    container:
        config["picard"]
    resources:
        java_mem_mb=1024,
        mem_mb=2000
    shell:
        """
        picard MarkDuplicates \
            -Xmx{resources.java_mem_mb}M \
            {params.extra} \
            INPUT={input[0]} \
            OUTPUT={output.bam} \
            METRICS_FILE={output.metrics} \
                2> {log}
        """

rule samtools_index:
    input:
        rules.mark_duplicates.output,
    output:
        os.path.join(
            config["workdir"],
            "bams/hdrr/bwamem2/marked/{run}/{sample}.bam.bai"
        ),
    log:
        os.path.join(
            config["workdir"], 
            "logs/samtools_index/hdrr/{run}/{sample}.log"
        ),
    resources:
        mem_mb = 100
    container:
        config["samtools"]
    shell:
        """
        samtools index \
            {input[0]} \
            {output[0]} \
                2> {log}
        """
