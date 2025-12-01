include: "common.smk"


rule plot_lionelmssq_prediction:
    input:
        pred_fragments="results/lionelmssq/{modus}/{seq}/{n_fragments}.tsv",
        pred_seq="results/lionelmssq/{modus}/{seq}/{n_fragments}.fasta",
        sim="data/{modus}/{seq}/{n_fragments}.tsv",
    output:
        "results/plots/lionelmssq_prediction/{modus}/{seq}/{n_fragments}.html",
    log:
        "logs/plots/lionelmssq_prediction/{modus}/{seq}/{n_fragments}.log",
    benchmark:
        "benchmarks/plots/lionelmssq_prediction/{modus}/{seq}/{n_fragments}.tsv"
    conda:
        "../envs/lionelmssq.yaml"
    threads: 1
    script:
        "../scripts/plot_lionelmssq_prediction.py"


rule evaluate_custom_simulation:
    input:
        collect_custom_simulations(
            "results/lionelmssq/simulation/{seq}/{n_fragments}.fasta"
        ),
    output:
        "results/evaluation/custom_simulation.tsv",
    log:
        "logs/evaluation/custom_simulation.log",
    benchmark:
        "benchmarks/evaluation/custom_simulation.tsv"
    conda:
        "../envs/lionelmssq.yaml"
    threads: 1
    script:
        "../scripts/evaluate_prediction.py"


rule plot_evaluation_custom_simulation:
    input:
        "results/evaluation/custom_simulation.tsv",
    output:
        donut="results/plots/evaluation/custom_simulation.donut.html",
        bar="results/plots/evaluation/custom_simulation.bar.html",
    log:
        "logs/plots/evaluation/custom_simulation.log",
    benchmark:
        "benchmarks/plots/evaluation/custom_simulation.tsv"
    conda:
        "../envs/lionelmssq.yaml"
    threads: 1
    script:
        "../scripts/plot_evaluation.py"


rule evaluate_random_simulation:
    input:
        collect_random_simulations(
            "results/lionelmssq/simulation/{seq}/{n_fragments}.fasta"
        ),
    output:
        "results/evaluation/random_simulation.tsv",
    log:
        "logs/evaluation/random_simulation.log",
    benchmark:
        "benchmarks/evaluation/random_simulation.tsv"
    conda:
        "../envs/lionelmssq.yaml"
    threads: 1
    script:
        "../scripts/evaluate_prediction.py"


rule plot_evaluation_random_simulation:
    input:
        "results/evaluation/random_simulation.tsv",
    output:
        donut="results/plots/evaluation/random_simulation.donut.html",
        bar="results/plots/evaluation/random_simulation.bar.html",
    log:
        "logs/plots/evaluation/random_simulation.log",
    benchmark:
        "benchmarks/plots/evaluation/random_simulation.tsv"
    conda:
        "../envs/lionelmssq.yaml"
    threads: 1
    script:
        "../scripts/plot_evaluation.py"


rule evaluate_experiment:
    input:
        collect_experiments("results/lionelmssq/experiment/{seq}/{n_fragments}.fasta"),
    output:
        "results/evaluation/experiment.tsv",
    log:
        "logs/evaluation/experiment.log",
    benchmark:
        "benchmarks/evaluation/experiment.tsv"
    conda:
        "../envs/lionelmssq.yaml"
    threads: 1
    script:
        "../scripts/evaluate_prediction.py"


rule plot_evaluation_experiment:
    input:
        "results/evaluation/experiment.tsv",
    output:
        donut="results/plots/evaluation/experiment.donut.html",
        bar="results/plots/evaluation/experiment.bar.html",
    log:
        "logs/plots/evaluation/experiment.log",
    benchmark:
        "benchmarks/plots/evaluation/experiment.tsv"
    conda:
        "../envs/lionelmssq.yaml"
    threads: 1
    script:
        "../scripts/plot_evaluation.py"
