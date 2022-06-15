rule merge_datasets:
    input:
        expand(os.path.join(
            config["rawdir"],
            "{{interval}}/{dataset}.csv"
            ),
            dataset = config["datasets"]
        ),
    output:
        os.path.join(
            config['workdir'],
            "merged/{interval}.csv"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/merge_datasets/{interval}.log"
        ),
    resources:
        mem_mb = 50000
    script:
        "../scripts/merge_datasets.py"