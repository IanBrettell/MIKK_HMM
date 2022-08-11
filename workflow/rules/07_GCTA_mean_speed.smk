## Create .phen files
## as specified here: https://gcta.freeforums.net/thread/247/greml-estimating-variance-explained-snps
rule create_phen_mean_speed:
    input:
        fam = expand(rules.create_bed_all.output.fam,
            bin_length = 5000,
            cov = 0.8
        ),
        phenos = os.path.join(
            config["workdir"],
            "mean_speed_F2/{interval}/{dge_sge}/{transformation}.csv"
            ),
        samples_file = config["F2_samples_file"]
    output:
        os.path.join(
            config["workdir"],
            "phens_ms/true/{interval}/{dge_sge}/{transformation}/{assay}.phen"
        ),
    log:
        os.path.join(
            config["workdir"],
            "logs/create_phen_mean_speed/{interval}/{dge_sge}/{transformation}/{assay}.log"
        ),
    params:
        assay = "{assay}"        
    resources:
        mem_mb = 20000,
    container:
        config["R_4.2.0"]
    script:
        "../scripts/create_phen_mean_speed.R"

rule permute_phen_mean_speed:
    input:
        fam = rules.create_phen_mean_speed.input.fam,
        phen = rules.create_phen_mean_speed.output,
    output:
        os.path.join(
            config["workdir"],
            "phens_ms/permuted/{interval}/{dge_sge}/{transformation}/{assay}/{seed}.phen"
        ),
    log:
        os.path.join(
            config["workdir"],
            "logs/permute_phen_mean_speed/{interval}/{dge_sge}/{transformation}/{assay}/{seed}.log"
        ),
    params:
        seed = "{seed}"
    resources:
        mem_mb = 500,
    container:
        config["R_4.2.0"]
    script:
        "../scripts/permute_phen.R"

# Set rule order for creating covariate files and running mlma-loco
# Because the covariate files aren't in the mlma rules' inputs
#ruleorder: create_covar > permute_covar > run_mlma_loco > run_mlma_loco_permuted

# Create covariate files
# as specified here: https://gcta.freeforums.net/thread/247/greml-estimating-variance-explained-snps
rule create_covar_mean_speed:
    input:
        fam = rules.create_phen_mean_speed.input.fam,
        samples_file = config["F2_samples_file"],
        covars_qc = config["covariate_cat_quant"]
    output:
        cat = os.path.join(
            config["workdir"],
            "covars_ms/true/{covars}.covar"
        ),
        quant = os.path.join(
            config["workdir"],
            "covars_ms/true/{covars}.qcovar"
        ),
    log:
        os.path.join(
            config["workdir"],
            "logs/create_covar_mean_speed/{covars}.log"
        ),
    params:
        covars = "{covars}"
    resources:
        mem_mb = 15000,
    container:
        config["R_4.2.0"]
    script:
        "../scripts/create_covar.R"

rule permute_covars_mean_speed:
    input:
        fam = rules.create_phen_mean_speed.input.fam,
        covar_cat = rules.create_covar_mean_speed.output.cat,
        covar_quant = rules.create_covar_mean_speed.output.quant,
    output:
        cat = os.path.join(
            config["workdir"],
            "covars_ms/permuted/{covars}/{seed}.covar"
        ),
        quant = os.path.join(
            config["workdir"],
            "covars_ms/permuted/{covars}/{seed}.qcovar"
        ),
    log:
        os.path.join(
            config["workdir"],
            "logs/permute_covars_mean_speed/{covars}/{seed}.log"
        ),
    params:
        covars = "{covars}",
        seed = "{seed}"
    resources:
        mem_mb = 15000,
    container:
        config["R_4.2.0"]
    script:
        "../scripts/permute_covars.R"

def set_covars(wildcards):
    if wildcards.covars == "None":
        out = ""
    else:
        covars_file_cat = os.path.join(
            config["workdir"], "covars_ms/true/" + wildcards.covars + ".covar")
        covars_file_quant = os.path.join(
            config["workdir"], "covars_ms/true/" + wildcards.covars + ".qcovar")
        out = '--covar ' + covars_file_cat + ' --qcovar ' + covars_file_quant
    return(out)

rule run_mlma_loco_mean_speed:
    input:
        bed = rules.create_bed_all.output.bed,
        phen = rules.create_phen_mean_speed.output,
        grm = rules.make_grm_loco_man.output,
    output:
        os.path.join(
            config["workdir"],
            "gcta/mlma_loco_mean_speed/true/hdrr/{interval}/{bin_length}/{cov}/{dge_sge}/{transformation}/{assay}/{covars}/{contig}.mlma"
        ),
    log:
        os.path.join(
            config["workdir"],
            "logs/run_mlma_loco/hdrr/{interval}/{bin_length}/{cov}/{dge_sge}/{transformation}/{assay}/{covars}/{contig}.log"
        ),
    params:
        in_pref = lambda wildcards, input: input.bed.replace(".bed", ""),
        grm_pref = lambda wildcards, input: input.grm[0].replace(".grm.bin", ""),
        out_pref = lambda wildcards, output: output[0].replace(".mlma", ""),
        contig = "{contig}",
        covars = set_covars,
    resources:
        mem_mb = 10000,
        threads = 1
    container:
        config["GCTA"]
    shell:
        """
        gcta64 \
            --mlma \
            --bfile {params.in_pref} \
            --grm {params.grm_pref} \
            --pheno {input.phen} \
            --chr {params.contig} \
            {params.covars} \
            --out {params.out_pref} \
            --autosome-num 24 \
            --thread-num {resources.threads} \
                2> {log}
        """

rule combine_mlma_outputs_mean_speed:
    input:
        expand(os.path.join(
            config["workdir"],
            "gcta/mlma_loco_mean_speed/true/hdrr/{{interval}}/{{bin_length}}/{{cov}}/{{dge_sge}}/{{transformation}}/{{assay}}/{{covars}}/{contig}.mlma"
            ),
                contig = list(range(1,22))
        ),
    output:
        os.path.join(
            config["workdir"],
            "gcta/mlma_loco_mean_speed_consol/true/hdrr/{interval}/{bin_length}/{cov}/{dge_sge}/{transformation}/{assay}/{covars}.mlma"
            ),
    log:
        os.path.join(
            config["workdir"],
            "logs/combine_mlma_outputs_mean_speed/hdrr/{interval}/{bin_length}/{cov}/{dge_sge}/{transformation}/{assay}/{covars}.log"
        ),
    resources:
        mem_mb = 10000,
    container:
        config["R_4.2.0"]
    script:
        "../scripts/combine_mlma_outputs.R"

def set_covars_permuted(wildcards):
    if wildcards.covars == "None":
        out = ""
    else:
        covars_file_cat = os.path.join(
            config["workdir"], 
            "covars_ms/permuted/" + wildcards.covars + "/" + wildcards.seed + ".covar"
            )
        covars_file_quant = os.path.join(
            config["workdir"], 
            "covars_ms/permuted/" + wildcards.covars + "/" + wildcards.seed + ".qcovar"
            )
        out = '--covar ' + covars_file_cat + ' --qcovar ' + covars_file_quant
    return(out)

rule run_mlma_loco_mean_speed_permuted:
    input:
        bed = rules.create_bed_all.output.bed,
        phen = rules.permute_phen_mean_speed.output,
        grm = rules.make_grm_loco_man.output,
    output:
        os.path.join(
            config["workdir"],
            "gcta/mlma_loco_mean_speed/permuted/hdrr/{interval}/{bin_length}/{cov}/{dge_sge}/{transformation}/{assay}/{covars}/{contig}/{seed}.mlma"
        ),
    log:
        os.path.join(
            config["workdir"],
            "logs/run_mlma_loco_mean_speed_permuted/hdrr/{interval}/{bin_length}/{cov}/{dge_sge}/{transformation}/{assay}/{covars}/{contig}/{seed}.log"
        ),
    params:
        in_pref = lambda wildcards, input: input.bed.replace(".bed", ""),
        grm_pref = lambda wildcards, input: input.grm[0].replace(".grm.bin", ""),
        out_pref = lambda wildcards, output: output[0].replace(".mlma", ""),
        covars = set_covars_permuted,
        contig = "{contig}"
    resources:
        mem_mb = 10000,
        threads = 1
    container:
        config["GCTA"]
    shell:
        """
        gcta64 \
            --mlma \
            --bfile {params.in_pref} \
            --grm {params.grm_pref} \
            --pheno {input.phen} \
            --chr {params.contig} \
            {params.covars} \
            --out {params.out_pref} \
            --autosome-num 24 \
            --thread-num {resources.threads} \
                2> {log}
        """

rule combine_mlma_outputs_mean_speed_permuted:
    input:
        expand(os.path.join(
            config["workdir"],
            "gcta/mlma_loco_mean_speed/permuted/hdrr/{{interval}}/{{bin_length}}/{{cov}}/{{dge_sge}}/{{transformation}}/{{assay}}/{{covars}}/{contig}/{{seed}}.mlma"
            ),
                contig = list(range(1,22))
        ),
    output:
        os.path.join(
            config["workdir"],
            "gcta/mlma_loco_mean_speed_consol/permuted/hdrr/{interval}/{bin_length}/{cov}/{dge_sge}/{transformation}/{assay}/{covars}/{seed}.mlma"
            ),
    log:
        os.path.join(
            config["workdir"],
            "logs/combine_mlma_outputs_mean_speed_permuted/hdrr/{interval}/{bin_length}/{cov}/{dge_sge}/{transformation}/{assay}/{covars}/{seed}.log"
        ),
    resources:
        mem_mb = 30000,
    container:
        config["R_4.2.0"]
    script:
        "../scripts/combine_mlma_outputs.R"

rule get_min_p_perms_mean_speed:
    input:
        expand(os.path.join(
            config["workdir"],
            "gcta/mlma_loco_mean_speed_consol/permuted/hdrr/{{interval}}/{{bin_length}}/{{cov}}/{{dge_sge}}/{{transformation}}/{{assay}}/{{covars}}/{seed}.mlma"
            ),
                seed = PERM_SEEDS         
        ),
    output:
        os.path.join(
            config["workdir"],
            "gcta/mlma_loco/min_p/hdrr/{interval}/{bin_length}/{cov}/{dge_sge}/{transformation}/{assay}/{covars}.csv"
        ),
    log:
        os.path.join(
            config["workdir"],
            "logs/get_min_p_perms_mean_speed/hdrr/{interval}/{bin_length}/{cov}/{dge_sge}/{transformation}/{assay}/{covars}.log"
        ),
    resources:
        mem_mb = 30000
    container:
        config["R_4.2.0"]
    script:
        "../scripts/get_min_p_perms.R"

rule get_manhattan_gcta_mean_speed:
    input:
        res = rules.combine_mlma_outputs_mean_speed.output,
        min_p = rules.get_min_p_perms_mean_speed.output,
    output:
        man = "book/figs/gcta/hdrr/mean_speed/{interval}/{bin_length}/{cov}/{dge_sge}/{transformation}/{assay}/{covars}.png",
        sig = "results/{interval}/{bin_length}/{cov}/{dge_sge}/{transformation}/{assay}/{covars}.tsv"
    log:
        os.path.join(
            config["workdir"],
            "logs/get_manhattan_gcta/hdrr/{interval}/{bin_length}/{cov}/{dge_sge}/{transformation}/{assay}/{covars}.log"
        ),    
    params:
        bin_length = "{bin_length}",
        cov = "{cov}",
        dge_sge = "{dge_sge}",
        transformation = "{transformation}",
        assay = "{assay}",
        covars = "{covars}",
        pheno = "mean_speed"
    resources:
        mem_mb = 30000
    container:
        config["R_4.1.3"]
    script:
        "../scripts/get_manhattan_gcta.R"
