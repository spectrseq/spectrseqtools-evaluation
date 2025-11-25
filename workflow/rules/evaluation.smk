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
    script:
        "../scripts/plot_lionelmssq_prediction.py"


rule plot_evaluation_custom_simulation:
    input:
        collect_custom_simulations(
            "results/lionelmssq/simulation/{seq}/{n_fragments}.fasta"
        ),
    output:
        donut="results/plots/evaluation/custom_simulation.donut.html",
        bar="results/plots/evaluation/custom_simulation.bar.html",
    log:
        "logs/plots/evaluation/custom_simulation.log",
    benchmark:
        "benchmarks/plots/evaluation/custom_simulation.tsv"
    conda:
        "../envs/lionelmssq.yaml"
    script:
        "../scripts/plot_evaluation.py"


rule plot_evaluation_random_simulation:
    input:
        collect_random_simulations(
            "results/lionelmssq/simulation/{seq}/{n_fragments}.fasta"
        ),
    output:
        donut="results/plots/evaluation/random_simulation.donut.html",
        bar="results/plots/evaluation/random_simulation.bar.html",
    log:
        "logs/plots/evaluation/random_simulation.log",
    benchmark:
        "benchmarks/plots/evaluation/random_simulation.tsv"
    conda:
        "../envs/lionelmssq.yaml"
    script:
        "../scripts/plot_evaluation.py"


rule plot_evaluation_experiment:
    input:
        collect_experiments("results/lionelmssq/experiment/{seq}/{n_fragments}.fasta"),
    output:
        donut="results/plots/evaluation/experiment.donut.html",
        bar="results/plots/evaluation/experiment.bar.html",
    log:
        "logs/plots/evaluation/experiment.log",
    benchmark:
        "benchmarks/plots/evaluation/experiment.tsv"
    conda:
        "../envs/lionelmssq.yaml"
    script:
        "../scripts/plot_evaluation.py"
