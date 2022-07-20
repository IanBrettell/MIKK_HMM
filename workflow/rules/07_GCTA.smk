# Create .ped and .map files with format
rule create_ped:
    input:
        expand(os.path.join(
            config["workdir"],
            "F2_with_genos/hdrr/{{bin_length}}/{{cov}}/{sample}_{pat}_{mat}.csv"
            ),
                zip,
                sample = SAMPLES_F2_ZIP,
                pat = PAT_ZIP,
                mat = MAT_ZIP
        ),
    output:
        ped = os.path.join(
            config["workdir"],
            "peds/F2/hdrr/{bin_length}/{cov}.ped"
        ),
        map = os.path.join(
            config["workdir"],
            "peds/F2/hdrr/{bin_length}/{cov}.map"
        ),
    log:
        os.path.join(
            config["workdir"],
            "logs/create_ped/{bin_length}/{cov}.log"
        ),
    resources:
        mem_mb = 30000,
    container:
        # requires tidyr >= v1.2
        config["tidyverse_4.1.3"]
    script:
        "../scripts/create_ped.R"

# Convert .ped to .bed
#rule create_bed:
#    input:
#        rules.create_ped.output.ped
#    output:
#        bed = os.path.join(
#            config["working_dir"],
#            "beds/F2/{ref}/{max_reads}/{bin_length}/{cov}.bed"
#        ),
#    log:
#        os.path.join(
#            config["working_dir"],
#            "logs/create_bed/{ref}/{max_reads}/{bin_length}/{cov}.log"
#        ),
#    params:
#        in_pref = lambda wildcards, input: input[0].replace(".ped", ""),
#        out_pref = lambda wildcards, output: output.bed.replace(".bed", ""),
#    resources:
#        mem_mb = 200
#    container:
#        config["plink1.9"]
#    shell:
#        """
#        plink1.9 \
#            --make-bed \
#            --no-fid \
#            --no-parents \
#            --no-sex \
#            --no-pheno \
#            --chr-set 24 no-xy no-mt \
#            --file {params.in_pref} \
#            --out {params.out_pref} \
#                2> {log}
#        """
#