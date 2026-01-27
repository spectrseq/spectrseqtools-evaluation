rule lionelmssq_experiment:
    input:
        fragments="data/experiment/{seq}/{n_fragments}.raw",
        meta="data/experiment/{seq}/{n_fragments}.meta.yaml",
    output:
        fragments="data/experiment/{seq}/{n_fragments}.tsv",
        singletons="data/experiment/{seq}/{n_fragments}.singletons.tsv",
        predictions="results/lionelmssq/experiment/{seq}/{n_fragments}.tsv",
        sequence="results/lionelmssq/experiment/{seq}/{n_fragments}.fasta",
    params:
        solver=config["solver"],
    log:
        "logs/lionelmssq/experiment/{seq}/{n_fragments}.log",
    benchmark:
        "benchmarks/lionelmssq/experiment/{seq}/{n_fragments}.tsv"
    conda:
        "../envs/lionelmssq.yaml"
    threads: 1
    shell:
        "lionelmssq --fragments {input.fragments} --meta {input.meta} "
        "--fragment-predictions {output.predictions} "
        "--sequence-prediction {output.sequence} "
        "--sequence-name 'lionelmssq_prediction_from_sim_{wildcards.seq}' "
        "--solver {params.solver} "
        "--threads {threads} "
        "2> {log}"


rule lionelmssq_simulation:
    input:
        fragments="data/simulation/{seq}/{n_fragments}.tsv",
        meta="data/simulation/{seq}/{n_fragments}.meta.yaml",
        singletons="data/simulation/{seq}/{n_fragments}.singletons.tsv",
    output:
        predictions="results/lionelmssq/simulation/{seq}/{n_fragments}.tsv",
        sequence="results/lionelmssq/simulation/{seq}/{n_fragments}.fasta",
    params:
        solver=config["solver"],
    log:
        "logs/lionelmssq/simulation/{seq}/{n_fragments}.log",
    benchmark:
        "benchmarks/lionelmssq/simulation/{seq}/{n_fragments}.tsv"
    conda:
        "../envs/lionelmssq.yaml"
    threads: 1
    shell:
        "lionelmssq --fragments {input.fragments} --meta {input.meta} "
        "--fragment-predictions {output.predictions} "
        "--sequence-prediction {output.sequence} "
        "--sequence-name 'lionelmssq_prediction_from_sim_{wildcards.seq}' "
        "--solver {params.solver} "
        "--threads {threads} "
        "2> {log}"


rule lionelmssq_comparison_study:
    input:
        fragments="comparison_study/{parameter}/{value}/{seq}/sample.tsv",
        meta="comparison_study/{parameter}/{value}/{seq}/sample.meta.yaml",
    output:
        predictions="results/comparison_study/{parameter}/{value}/{seq}/sample.tsv",
        sequence="results/comparison_study/{parameter}/{value}/{seq}/sample.fasta",
    params:
        solver=config["solver"],
    log:
        "logs/comparison_study/{parameter}/{value}/{seq}/sample.log",
    benchmark:
        "benchmarks/comparison_study/{parameter}/{value}/{seq}/sample.tsv"
    conda:
        "../envs/lionelmssq.yaml"
    threads: 1
    shell:
        "lionelmssq --fragments {input.fragments} --meta {input.meta} "
        "--fragment-predictions {output.predictions} "
        "--sequence-prediction {output.sequence} "
        "--sequence-name 'lionelmssq_prediction_from_sim_{wildcards.seq}' "
        "--solver {params.solver} "
        "--threads {threads} "
        "2> {log}"
