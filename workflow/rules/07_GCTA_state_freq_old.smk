## Create .phen files
## as specified here: https://gcta.freeforums.net/thread/247/greml-estimating-variance-explained-snps
#rule create_phen:
#    input:
#        ped = expand(rules.create_ped.output.ped,
#            bin_length = 5000,
#            cov = 0.8
#        ),
#        phenos = os.path.join(
#            config["workdir"],
#            "state_freq_F2/dist_angle/0.05/15/{dge_sge}/{transformation}/{state}.csv"
#            ),
#        samples_file = config["F2_samples_file"]
#    output:
#        os.path.join(
#            config["workdir"],
#            "phens/true/{dge_sge}/{transformation}/{assay}/{state}.phen"
#        ),
#    log:
#        os.path.join(
#            config["workdir"],
#            "logs/create_phen/{dge_sge}/{transformation}/{assay}/{state}.log"
#        ),
#    params:
#        assay = "{assay}"        
#    resources:
#        mem_mb = 20000,
#    container:
#        # requires tidyr >= v1.2
#        config["tidyverse_4.1.3"]
#    script:
#        "../scripts/create_phen.R"
#
#rule permute_phen:
#    input:
#        rules.create_phen.output,
#    output:
#        os.path.join(
#            config["workdir"],
#            "phens/permuted/{dge_sge}/{transformation}/{assay}/{state}/{seed}.phen"
#        ),
#    log:
#        os.path.join(
#            config["workdir"],
#            "logs/permute_phen/{dge_sge}/{transformation}/{assay}/{state}/{seed}.log"
#        ),
#    params:
#        phenotype = "state_freq",
#        seed = "{seed}"
#    resources:
#        mem_mb = 500,
#    container:
#        # requires tidyr >= v1.2
#        config["tidyverse_4.1.3"]
#    script:
#        "../scripts/permute_phen.R"
#
## Set rule order for creating covariate files and running mlma-loco
## Because the covariate files aren't in the mlma rules' inputs
#ruleorder: create_covar > run_mlma_loco > run_mlma_loco_permuted
#
## Create covariate files
## as specified here: https://gcta.freeforums.net/thread/247/greml-estimating-variance-explained-snps
#rule create_covar:
#    input:
#        #genos = rules.process_rc_blocks.output,
#        ped = expand(rules.create_ped.output.ped,
#            bin_length = 5000,
#            cov = 0.8
#        ),
#        samples_file = config["F2_samples_file"]
#    output:
#        cat = os.path.join(
#            config["workdir"],
#            "covars/true/{covars}.covar"
#        ),
#        quant = os.path.join(
#            config["workdir"],
#            "covars/true/{covars}.qcovar"
#        ),
#    log:
#        os.path.join(
#            config["workdir"],
#            "logs/create_covar/{covars}.log"
#        ),
#    params:
#        covars = "{covars}"
#    resources:
#        mem_mb = 15000,
#    container:
#        # requires tidyr >= v1.2
#        config["tidyverse_4.1.3"]
#    script:
#        "../scripts/create_covar.R"
#
#rule permute_covars:
#    input:
#        ped = expand(rules.create_ped.output.ped,
#            bin_length = 5000,
#            cov = 0.8
#        ),
#        covar_cat = os.path.join(
#            config["workdir"],
#            "covars/true/All.covar"
#        ),
#        covar_quant = os.path.join(
#            config["workdir"],
#            "covars/true/All.qcovar"
#        ),
#    output:
#        cat = os.path.join(
#            config["workdir"],
#            "covars/permuted/{covars}/{seed}.covar"
#        ),
#        quant = os.path.join(
#            config["workdir"],
#            "covars/permuted/{covars}/{seed}.qcovar"
#        ),
#    log:
#        os.path.join(
#            config["workdir"],
#            "logs/permute_covar/{covars}/{seed}.log"
#        ),
#    params:
#        covars = "{covars}",
#        seed = "{seed}"
#    resources:
#        mem_mb = 15000,
#    container:
#        config["R_4.2.0"]
#    script:
#        "../scripts/permute_covars.R"
#
#def set_covars(wildcards):
#    if wildcards.covars == "None":
#        out = ""
#    else:
#        covars_file_cat = os.path.join(
#            config["workdir"], "covars/true/" + wildcards.covars + ".covar")
#        covars_file_quant = os.path.join(
#            config["workdir"], "covars/true/" + wildcards.covars + ".qcovar")
#        out = '--covar ' + covars_file_cat + ' --qcovar ' + covars_file_quant
#    return(out)
#
#rule run_mlma_loco:
#    input:
#        bed = rules.create_bed.output.bed,
#        phen = rules.create_phen.output,
#        #excl_samples = rules.create_excluded_samples_list.output,
#    output:
#        os.path.join(
#            config["workdir"],
#            "gcta/mlma_loco/true/hdrr/{bin_length}/{cov}/{dge_sge}/{transformation}/{assay}/{state}/{covars}.loco.mlma"
#        ),
#    log:
#        os.path.join(
#            config["workdir"],
#            "logs/run_mlma_loco/hdrr/{bin_length}/{cov}/{dge_sge}/{transformation}/{assay}/{state}/{covars}.log"
#        ),
#    params:
#        in_pref = lambda wildcards, input: input.bed.replace(".bed", ""),
#        out_pref = lambda wildcards, output: output[0].replace(".loco.mlma", ""),
#        covars = set_covars,
#    resources:
#        mem_mb = 5000,
#        threads = 1
#    container:
#        config["GCTA"]
#    shell:
#        """
#        gcta64 \
#            --mlma-loco \
#            --bfile {params.in_pref} \
#            --pheno {input.phen} \
#            {params.covars} \
#            --out {params.out_pref} \
#            --autosome-num 24 \
#            --thread-num {resources.threads} \
#                2> {log}
#        """
#
#def set_covars_permuted(wildcards):
#    if wildcards.covars == "None":
#        out = ""
#    else:
#        covars_file_cat = os.path.join(
#            config["workdir"], 
#            "covars/permuted/" + wildcards.covars + "/" + wildcards.seed + ".covar"
#            )
#        covars_file_quant = os.path.join(
#            config["workdir"], 
#            "covars/permuted/" + wildcards.covars + "/" + wildcards.seed + ".qcovar"
#            )
#        out = '--covar ' + covars_file_cat + ' --qcovar ' + covars_file_quant
#    return(out)
#
#rule run_mlma_loco_permuted:
#    input:
#        bed = rules.create_bed.output.bed,
#        phen = rules.permute_phen.output,
#    output:
#        os.path.join(
#            config["workdir"],
#            "gcta/mlma_loco/permuted/hdrr/{bin_length}/{cov}/{dge_sge}/{transformation}/{assay}/{state}/{covars}/{seed}.loco.mlma"
#        ),
#    log:
#        os.path.join(
#            config["workdir"],
#            "logs/run_mlma_loco_permuted/hdrr/{bin_length}/{cov}/{dge_sge}/{transformation}/{assay}/{state}/{covars}/{seed}.log"
#        ),
#    params:
#        in_pref = lambda wildcards, input: input.bed.replace(".bed", ""),
#        out_pref = lambda wildcards, output: output[0].replace(".loco.mlma", ""),
#        covars = set_covars_permuted,
#    resources:
#        mem_mb = 5000,
#        threads = 1
#    container:
#        config["GCTA"]
#    shell:
#        """
#        gcta64 \
#            --mlma-loco \
#            --bfile {params.in_pref} \
#            --pheno {input.phen} \
#            {params.covars} \
#            --out {params.out_pref} \
#            --autosome-num 24 \
#            --thread-num {resources.threads} \
#                2> {log}
#        """
#
#rule get_min_p_perms:
#    input:
#        expand(os.path.join(
#            config["workdir"],
#            "gcta/mlma_loco/permuted/hdrr/{{bin_length}}/{{cov}}/{{dge_sge}}/{{transformation}}/{{assay}}/{{state}}/{{covars}}/{seed}.loco.mlma"
#            ),
#                seed = PERM_SEEDS         
#        ),
#    output:
#        os.path.join(
#            config["workdir"],
#            "gcta/mlma_loco/min_p/hdrr/{bin_length}/{cov}/{dge_sge}/{transformation}/{assay}/{state}/{covars}.csv"
#        ),
#    log:
#        os.path.join(
#            config["workdir"],
#            "logs/get_min_p_perms/hdrr/{bin_length}/{cov}/{dge_sge}/{transformation}/{assay}/{state}/{covars}.log"
#        ),
#    resources:
#        mem_mb = 8000
#    container:
#        config["R_4.2.0"]
#    script:
#        "../scripts/get_min_p_perms.R"
#
#rule get_manhattan_gcta:
#    input:
#        res = rules.run_mlma_loco.output,
#        min_p = rules.get_min_p_perms.output,
#    output:
#        man = "book/figs/gcta/hdrr/{bin_length}/{cov}/{dge_sge}/{transformation}/{assay}/{state}_{covars}.png"
#    log:
#        os.path.join(
#            config["workdir"],
#            "logs/get_manhattan_gcta/hdrr/{bin_length}/{cov}/{dge_sge}/{transformation}/{assay}/{state}_{covars}.log"
#        ),    
#    params:
#        bin_length = "{bin_length}",
#        cov = "{cov}",
#        dge_sge = "{dge_sge}",
#        transformation = "{transformation}",
#        assay = "{assay}",
#        state = "{state}",
#        covars = "{covars}"
#    resources:
#        mem_mb = 6000
#    container:
#        config["R_4.1.3"]
#    script:
#        "../scripts/get_manhattan_gcta.R"
#
