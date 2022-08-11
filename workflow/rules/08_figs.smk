#rule karyoplots:
#    input:
#        data = rules.true_hmmlearn.output.csv,
#        line_cols = rules.get_line_ranks_and_colours.output.csv
#    output:
#        karyoplot_no_missing = "book/figs/hmmlearn_true/{bin_length}/{cov}/karyoplot_no_missing.png",
#        karyoplot_with_missing = "book/figs/hmmlearn_true/{bin_length}/{cov}/karyoplot_wi_missing.png",
#    log:
#        os.path.join(
#            config["workdir"],
#            "logs/karyoplots/{bin_length}/{cov}.log"
#        ),
#    params:
#        cov = "{cov}",
#        bin_length = "{bin_length}",
#    resources:
#        mem_mb = 50000
#    container:
#        config["R_4.1.3"]
#    script:
#        "../scripts/karyoplots.R"

rule compile_sig_mans:
    input:
        expand("book/figs/gcta/hdrr/state_freq/{{variables}}/{{interval}}/{{n_states}}/{{bin_length}}/{{cov}}/{{dge_sge}}/{{transformation}}/{{assay}}/{state}_{{covars}}.png",
            state = list(range(1,16)),
        ),
    output:
        pdf = "book/figs/mans_compiled/hdrr/{variables}/{interval}/{n_states}/{bin_length}/{cov}/{dge_sge}/{transformation}/{assay}/{dge_sge}_{assay}_covars-{covars}.pdf",
    log:
        os.path.join(
            config["workdir"],
            "logs/compile_sig_mans/{variables}/{interval}/{n_states}/{bin_length}/{cov}/{dge_sge}/{transformation}/{assay}/{covars}.log"
        ),
    params:
        sheet_id = config["aov_google_sheet_select"],
        dge_sge = "{dge_sge}"
    resources:
        mem_mb = 10000
    container:
        config["ImageMagick_7.1.0.43"]
    shell:
        """
        convert {input} {output.pdf} \
            2> {log}
        """ 

##################
# GRMs
##################

#rule plot_grm:
#    input:
#        grm = rules.make_grm.output,
#        F2_samples = config["F2_samples_file"]
#    output:
#        png = "book/figs/grm/{bin_length}/{cov}/grm.png"
#    log:
#        os.path.join(
#            config["workdir"],
#            "logs/plot_grm/{bin_length}/{cov}.log"
#        ),
#    params:
#        grm_pref = lambda wildcards, input: input.grm[0].replace(".grm.bin", ""),
#    resources:
#        mem_mb = 5000,
#    container:
#        config["R_4.2.0"]
#    script:
#        "../scripts/plot_grm.R"
#
#rule plot_grm_per_chrom:
#    input:
#        grm = rules.make_grm_per_chrom.output,
#        F2_samples = config["F2_samples_file"]
#    output:
#        png = "book/figs/grm_per_chrom/{bin_length}/{cov}/{contig}_grm.png"
#    log:
#        os.path.join(
#            config["workdir"],
#            "logs/plot_grm_per_chrom/{bin_length}/{cov}/{contig}.log"
#        ),
#    params:
#        grm_pref = lambda wildcards, input: input.grm[0].replace(".grm.bin", ""),
#    resources:
#        mem_mb = 5000,
#    container:
#        config["R_4.2.0"]
#    script:
#        "../scripts/plot_grm.R"
#
#rule plot_grm_inbred_per_chrom:
#    input:
#        grm = rules.make_grm_inbred_per_chrom.output,
#        F2_samples = config["F2_samples_file"]
#    output:
#        png = "book/figs/grm_inbred_per_chrom/{bin_length}/{cov}/{contig}_grm.png"
#    log:
#        os.path.join(
#            config["workdir"],
#            "logs/plot_grm_inbred_per_chrom/{bin_length}/{cov}/{contig}.log"
#        ),
#    params:
#        grm_pref = lambda wildcards, input: input.grm[0].replace(".grm.bin", ""),
#    resources:
#        mem_mb = 5000,
#    container:
#        config["R_4.2.0"]
#    script:
#        "../scripts/plot_grm.R"
#
#rule plot_grm_no_miss:
#    input:
#        grm = rules.make_grm_no_miss.output,
#        F2_samples = config["F2_samples_file"]
#    output:
#        png = "book/figs/grm_no_miss/{bin_length}/{cov}.png"
#    log:
#        os.path.join(
#            config["workdir"],
#            "logs/plot_grm_no_miss/{bin_length}/{cov}.log"
#        ),
#    params:
#        grm_pref = lambda wildcards, input: input.grm[0].replace(".grm.bin", ""),
#    resources:
#        mem_mb = 5000,
#    container:
#        config["R_4.2.0"]
#    script:
#        "../scripts/plot_grm.R"
#
#rule plot_grm_inbred_no_miss:
#    input:
#        grm = rules.make_grm_inbred_no_miss.output,
#        F2_samples = config["F2_samples_file"]
#    output:
#        png = "book/figs/grm_inbred_no_miss/{bin_length}/{cov}.png"
#    log:
#        os.path.join(
#            config["workdir"],
#            "logs/plot_grm_inbred_no_miss/{bin_length}/{cov}.log"
#        ),
#    params:
#        grm_pref = lambda wildcards, input: input.grm[0].replace(".grm.bin", ""),
#    resources:
#        mem_mb = 5000,
#    container:
#        config["R_4.2.0"]
#    script:
#        "../scripts/plot_grm.R"
#