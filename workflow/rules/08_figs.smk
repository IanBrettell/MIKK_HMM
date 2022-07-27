rule karyoplots:
    input:
        data = rules.true_hmmlearn.output.csv,
        line_cols = rules.get_line_ranks_and_colours.output.csv
    output:
        karyoplot_no_missing = "book/figs/hmmlearn_true/{bin_length}/{cov}/karyoplot_no_missing.png",
        karyoplot_with_missing = "book/figs/hmmlearn_true/{bin_length}/{cov}/karyoplot_wi_missing.png",
    log:
        os.path.join(
            config["workdir"],
            "logs/karyoplots/{bin_length}/{cov}.log"
        ),
    params:
        cov = "{cov}",
        bin_length = "{bin_length}",
    resources:
        mem_mb = 50000
    container:
        config["R_4.1.3"]
    script:
        "../scripts/karyoplots.R"

rule compile_sig_mans:
    input:
        expand("book/figs/gcta/hdrr/{{bin_length}}/{{cov}}/{{dge_sge}}/{{transformation}}/{{assay}}/{state}_{{covars}}.png",
            state = list(range(1,16)),
        ),
    output:
        pdf = "book/figs/mans_compiled/hdrr/{bin_length}/{cov}/{transformation}/{dge_sge}/{assay}/{dge_sge}_{assay}_covars-{covars}.pdf",
    log:
        os.path.join(
            config["workdir"],
            "logs/compile_sig_mans/{bin_length}/{cov}/{dge_sge}/{transformation}/{assay}/{covars}.log"
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

rule plot_grm:
    input:
        grm = rules.make_grm.output,
        grm_inbred = rules.make_grm_inbred.output,
    output:
        grm_png = "book/figs/grm/{bin_length}/{cov}/grm.png",
        grm_inbred_png = "book/figs/grm/{bin_length}/{cov}/grm_inbred.png",
    log:
        os.path.join(
            config["workdir"],
            "logs/plot_grm/{bin_length}/{cov}.log"
        ),
    params:
        grm_pref = lambda wildcards, input: input.grm.replace(".grm.bin", ""),
        grm_inbred_pref = lambda wildcards, input: input.grm_inbred.replace(".grm.bin", ""),
    resources:
        mem_mb = 5000,
    container:
        config["R_4.2.0"]
    script:
        "../scripts/plot_grm.R"
    
        