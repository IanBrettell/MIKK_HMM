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

# Get list of F2 samples and pat and mat lines (to zip)
df_f2 = pd.read_csv(config["F2_samples_file"])
SAMPLES_F2_ZIP = df_f2['finclip_id']
PAT_ZIP = df_f2['pat_line'].tolist()
MAT_ZIP = df_f2['mat_line'].tolist()

# Create new lists of unique combinations of paternal and maternal lines
## Combine lists
PAT_MAT = zip(PAT_ZIP, MAT_ZIP)
## Get unique combinations
PM_LIST = set(list(PAT_MAT))
## Write to new lists
PAT_UQ = [i[0] for i in PM_LIST]
MAT_UQ = [i[1] for i in PM_LIST]

# Get separate lists for F2 and Kiyosu closed-capture
F2_mask = [not i.startswith('K') for i in SAMPLES]
KCC_mask = [i.startswith('K') for i in SAMPLES]

SAMPLES_KCC = [i for i in SAMPLES if i.startswith('K')]

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
        ref = os.path.join(
            config["workdir"],
            "refs/hdrr.fasta"
        ),
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
        extra=lambda wildcards: "VALIDATION_STRINGENCY=LENIENT TMP_DIR=" + config["tmpdir"],
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
        extra = lambda wildcards: "REMOVE_DUPLICATES=true TMP_DIR=" + config["tmpdir"]
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

rule merge_bams:
    input:
        expand(os.path.join(
            config["workdir"],
            "bams/hdrr/bwamem2/marked/{run}/{{sample}}.bam"),
                run = config["runs"]
        ),
    output:
        bam = os.path.join(
            config["workdir"],
            "bams/hdrr/bwamem2/merged/{sample}.bam"
        ),
    log:
        os.path.join(
            config["workdir"],
            "logs/merge_bams/{sample}.log"
        ),
    params:
        extra = lambda wildcards: "VALIDATION_STRINGENCY=LENIENT TMP_DIR=" + config["tmpdir"],
        in_files = lambda wildcards, input: " I=".join(input)
    resources:
        java_mem_mb=1024,
        mem_mb = 5000
    container:
        config["picard"]
    shell:
        """
        picard MergeSamFiles \
            -Xmx{resources.java_mem_mb}M \
            {params.extra} \
            INPUT={params.in_files} \
            OUTPUT={output} \
                &> {log}
        """

rule samtools_index:
    input:
        rules.merge_bams.output,
    output:
        os.path.join(
            config["workdir"],
            "bams/hdrr/bwamem2/merged/{sample}.bam.bai"
        ),
    log:
        os.path.join(
            config["workdir"], 
            "logs/samtools_index/hdrr/{sample}.log"
        ),
    resources:
        mem_mb = 2000
    container:
        config["samtools"]
    shell:
        """
        samtools index \
            {input[0]} \
            {output[0]} \
                2> {log}
        """

rule get_coverage:
    input:
        bam = rules.merge_bams.output,
        ind = rules.samtools_index.output,
    output:
        os.path.join(
            config["workdir"],
            "coverage/hdrr/bwamem2/{sample}.txt"
        ),
    log:
        os.path.join(
            config["workdir"], 
            "logs/get_coverage/hdrr/{sample}.log"
        ),
    resources:
        mem_mb = 2000
    container:
        config["samtools"]
    shell:
        """
        samtools coverage \
            {input.bam} >\
                {output[0]} \
                    2> {log}
        """

rule plot_coverage:
    input:
        expand(rules.get_coverage.output,
            sample = SAMPLES
        )
    output:
        png = "book/figs/coverage/all.png",
        pdf = "book/figs/coverage/all.pdf"
    log:
        os.path.join(
            config["workdir"], 
            "logs/plot_coverage/all.log"
        ),
    resources:
        mem_mb = 2000
    container:
        config["R_4.2.0"]
    script:
        "../scripts/plot_coverage.R"

