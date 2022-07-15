# Reahead mikk vcf to convert cram IDs to line IDs
rule rehead_mikk_vcf:
    input:
        vcf = config["mikk_vcf"],
        rehead_file = "config/mikk_vcf_cram2line-id_dupes-ed.txt"
    output:
        os.path.join(
            config["workdir"],
            "mikk_vcf/mikk_rehead.vcf.gz"
        ),
    log:
        os.path.join(
            config["workdir"],
            "logs/rehead_mikk_vcf/all.log"
        )
    resources:
        mem_mb = 2000
    container:
        config["bcftools_1.14"]
    shell:
        """
        bcftools reheader \
            --samples {input.rehead_file} \
            --output {output[0]} \
            {input.vcf} \
                2> {log}
        """

# Extract genotypes of parental lines for all biallelic SNPs 
rule extract_parental_snps:
    input:
        rules.rehead_mikk_vcf.output,
    output:
        os.path.join(
            config["workdir"],
            "genos/F0/biallelic_snps/{pat}_{mat}.txt"
        ),
    log:
        os.path.join(
            config["workdir"],
            "logs/extract_parental_snps/{pat}_{mat}.log"
        ),
    params:
        pat = "{pat}",
        mat = "{mat}"
    resources:
        mem_mb = 2000
    container:
        config["bcftools_1.14"]
    shell:
        """
        bcftools view \
            --samples {params.pat},{params.mat} \
            --min-alleles 2 \
            --max-alleles 2 \
            --types snps \
            --output-type u \
            {input} |\
        bcftools query \
            --print-header \
            --format '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT\\t%AD]\\n' \
            --output {output} \
                2> {log}
        """

# Extract target SNPs for F2 mapping
## We want SNPs that are homozygous-divergent between the paternal and maternal lines
rule extract_homo_div_snps:
    input:
        genos = rules.extract_parental_snps.output,
    output:
        full = os.path.join(
            config["workdir"],
            "genos/F0/homo_divergent/{pat}_{mat}.csv"
        ),
        sites = os.path.join(
            config["workdir"],
            "sites_files/F0/homo_divergent/{pat}_{mat}.txt"
        ),
    log:
        os.path.join(
            config["workdir"],
            "logs/extract_homo_div_snps/{pat}_{mat}.log"
        ),
    params:
        # set minimum allele depth
        min_ad = 5,
        pat = "{pat}",
        mat = "{mat}"        
    container:
        config["R_4.2.0"]
    resources:
        mem_mb = 20000
    script:
        "../scripts/extract_homo_div_snps.R"
