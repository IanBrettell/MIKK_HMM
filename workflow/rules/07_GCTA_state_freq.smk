# Create .phen files
# as specified here: https://gcta.freeforusf.net/thread/247/greml-estimating-variance-explained-snps
rule create_phen_state_freq:
    input:
        fam = expand(rules.create_bed_all.output.fam,
            bin_length = 5000,
            cov = 0.8
        ),
        phenos = os.path.join(
            config["workdir"],
            "state_freq_F2/{variables}/{interval}/{n_states}/{dge_sge}/{transformation}/{state}.csv"
            ),
        samples_file = config["F2_samples_file"]
    output:
        os.path.join(
            config["workdir"],
            "phens_sf/true/{variables}/{interval}/{n_states}/{dge_sge}/{transformation}/{assay}/{state}.phen"
        ),
    log:
        os.path.join(
            config["workdir"],
            "logs/create_phen_state_freq/{variables}/{interval}/{n_states}/{dge_sge}/{transformation}/{assay}/{state}.log"
        ),
    params:
        assay = "{assay}"        
    resources:
        mem_mb = 20000,
    container:
        config["R_4.2.0"]
    script:
        "../scripts/create_phen_state_freq.R"

rule permute_phen_state_freq:
    input:
        fam = rules.create_phen_state_freq.input.fam,
        phen = rules.create_phen_state_freq.output,
    output:
        os.path.join(
            config["workdir"],
            "phens_sf/permuted/{variables}/{interval}/{n_states}/{dge_sge}/{transformation}/{assay}/{state}/{seed}.phen"
        ),
    log:
        os.path.join(
            config["workdir"],
            "logs/permute_phen_state_freq/{variables}/{interval}/{n_states}/{dge_sge}/{transformation}/{assay}/{state}/{seed}.log"
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
rule create_covar_state_freq:
    input:
        fam = rules.create_phen_state_freq.input.fam,
        samples_file = config["F2_samples_file"],
        covars_qc = config["covariate_cat_quant"]
    output:
        cat = os.path.join(
            config["workdir"],
            "covars_sf/true/{covars}.covar"
        ),
        quant = os.path.join(
            config["workdir"],
            "covars_sf/true/{covars}.qcovar"
        ),
    log:
        os.path.join(
            config["workdir"],
            "logs/create_covar_state_freq/{covars}.log"
        ),
    params:
        covars = "{covars}"
    resources:
        mem_mb = 15000,
    container:
        config["R_4.2.0"]
    script:
        "../scripts/create_covar.R"

rule permute_covars_state_freq:
    input:
        fam = rules.create_phen_state_freq.input.fam,
        covar_cat = rules.create_covar_state_freq.output.cat,
        covar_quant = rules.create_covar_state_freq.output.quant,
    output:
        cat = os.path.join(
            config["workdir"],
            "covars_sf/permuted/{covars}/{seed}.covar"
        ),
        quant = os.path.join(
            config["workdir"],
            "covars_sf/permuted/{covars}/{seed}.qcovar"
        ),
    log:
        os.path.join(
            config["workdir"],
            "logs/permute_covars_state_freq/{covars}/{seed}.log"
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
            config["workdir"], "covars_sf/true/" + wildcards.covars + ".covar")
        covars_file_quant = os.path.join(
            config["workdir"], "covars_sf/true/" + wildcards.covars + ".qcovar")
        out = '--covar ' + covars_file_cat + ' --qcovar ' + covars_file_quant
    return(out)

rule run_mlma_loco_state_freq:
    input:
        bed = rules.create_bed_all.output.bed,
        phen = rules.create_phen_state_freq.output,
        grm = rules.make_grm_loco_man.output,
    output:
        os.path.join(
            config["workdir"],
            "gcta/mlma_loco_state_freq/true/hdrr/{variables}/{interval}/{n_states}/{bin_length}/{cov}/{dge_sge}/{transformation}/{assay}/{state}/{covars}/{contig}.mlma"
        ),
    log:
        os.path.join(
            config["workdir"],
            "logs/run_mlma_loco/hdrr/{variables}/{interval}/{n_states}/{bin_length}/{cov}/{dge_sge}/{transformation}/{assay}/{state}/{covars}/{contig}.log"
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

rule combine_mlma_outputs_state_freq:
    input:
        expand(os.path.join(
            config["workdir"],
            "gcta/mlma_loco_state_freq/true/hdrr/{{variables}}/{{interval}}/{{n_states}}/{{bin_length}}/{{cov}}/{{dge_sge}}/{{transformation}}/{{assay}}/{{state}}/{{covars}}/{contig}.mlma"
            ),
                contig = list(range(1,22))
        ),
    output:
        os.path.join(
            config["workdir"],
            "gcta/mlma_loco_state_freq_consol/true/hdrr/{variables}/{interval}/{n_states}/{bin_length}/{cov}/{dge_sge}/{transformation}/{assay}/{state}/{covars}.mlma"
            ),
    log:
        os.path.join(
            config["workdir"],
            "logs/combine_mlma_outputs_state_freq/hdrr/{variables}/{interval}/{n_states}/{bin_length}/{cov}/{dge_sge}/{transformation}/{assay}/{state}/{covars}.log"
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
            "covars_sf/permuted/" + wildcards.covars + "/" + wildcards.seed + ".covar"
            )
        covars_file_quant = os.path.join(
            config["workdir"], 
            "covars_sf/permuted/" + wildcards.covars + "/" + wildcards.seed + ".qcovar"
            )
        out = '--covar ' + covars_file_cat + ' --qcovar ' + covars_file_quant
    return(out)

rule run_mlma_loco_state_freq_permuted:
    input:
        bed = rules.create_bed_all.output.bed,
        phen = rules.permute_phen_state_freq.output,
        grm = rules.make_grm_loco_man.output,
    output:
        os.path.join(
            config["workdir"],
            "gcta/mlma_loco_state_freq/permuted/hdrr/{variables}/{interval}/{n_states}/{bin_length}/{cov}/{dge_sge}/{transformation}/{assay}/{state}/{covars}/{contig}/{seed}.mlma"
        ),
    log:
        os.path.join(
            config["workdir"],
            "logs/run_mlma_loco_state_freq_permuted/hdrr/{variables}/{interval}/{n_states}/{bin_length}/{cov}/{dge_sge}/{transformation}/{assay}/{state}/{covars}/{contig}/{seed}.log"
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

rule combine_mlma_outputs_state_freq_permuted:
    input:
        expand(os.path.join(
            config["workdir"],
            "gcta/mlma_loco_state_freq/permuted/hdrr/{{variables}}/{{interval}}/{{n_states}}/{{bin_length}}/{{cov}}/{{dge_sge}}/{{transformation}}/{{assay}}/{{state}}/{{covars}}/{contig}/{{seed}}.mlma"
            ),
                contig = list(range(1,22))
        ),
    output:
        os.path.join(
            config["workdir"],
            "gcta/mlma_loco_state_freq_consol/permuted/hdrr/{variables}/{interval}/{n_states}/{bin_length}/{cov}/{dge_sge}/{transformation}/{assay}/{state}/{covars}/{seed}.mlma"
            ),
    log:
        os.path.join(
            config["workdir"],
            "logs/combine_mlma_outputs_state_freq_permuted/hdrr/{variables}/{interval}/{n_states}/{bin_length}/{cov}/{dge_sge}/{transformation}/{assay}/{state}/{covars}/{seed}.log"
        ),
    resources:
        mem_mb = 30000,
    container:
        config["R_4.2.0"]
    script:
        "../scripts/combine_mlma_outputs.R"

rule get_min_p_perms_state_freq:
    input:
        expand(os.path.join(
            config["workdir"],
            "gcta/mlma_loco_state_freq_consol/permuted/hdrr/{{variables}}/{{interval}}/{{n_states}}/{{bin_length}}/{{cov}}/{{dge_sge}}/{{transformation}}/{{assay}}/{{state}}/{{covars}}/{seed}.mlma"
            ),
                seed = PERM_SEEDS         
        ),
    output:
        os.path.join(
            config["workdir"],
            "gcta/mlma_loco/min_p/hdrr/{variables}/{interval}/{n_states}/{bin_length}/{cov}/{dge_sge}/{transformation}/{assay}/{state}/{covars}.csv"
        ),
    log:
        os.path.join(
            config["workdir"],
            "logs/get_min_p_perms_state_freq/hdrr/{variables}/{interval}/{n_states}/{bin_length}/{cov}/{dge_sge}/{transformation}/{assay}/{state}/{covars}.log"
        ),
    resources:
        mem_mb = 30000
    container:
        config["R_4.2.0"]
    script:
        "../scripts/get_min_p_perms.R"

rule get_manhattan_gcta_state_freq:
    input:
        res = rules.combine_mlma_outputs_state_freq.output,
        min_p = rules.get_min_p_perms_state_freq.output,
    output:
        man = "book/figs/gcta/hdrr/state_freq/{variables}/{interval}/{n_states}/{bin_length}/{cov}/{dge_sge}/{transformation}/{assay}/{state}_{covars}.png"
    log:
        os.path.join(
            config["workdir"],
            "logs/get_manhattan_gcta/hdrr/{variables}/{interval}/{n_states}/{bin_length}/{cov}/{dge_sge}/{transformation}/{assay}/{state}_{covars}.log"
        ),    
    params:
        bin_length = "{bin_length}",
        cov = "{cov}",
        dge_sge = "{dge_sge}",
        transformation = "{transformation}",
        assay = "{assay}",
        covars = "{covars}",
        pheno = "state_freq",
        state = "{state}"
    resources:
        mem_mb = 30000
    container:
        config["R_4.1.3"]
    script:
        "../scripts/get_manhattan_gcta.R"
