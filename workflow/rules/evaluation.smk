# TODO implement script
rule plot_lionelmssq_prediction:
    input:
        tsv="results/lionelmssq/{seq}/{n_fragments}.tsv",
        fasta="results/lionelmssq/{seq}/{n_fragments}.fasta",
    output:
        "results/plots/lionelmssq_prediction/{seq}/{n_fragments}.html",
    log:
        "logs/plots/lionelmssq_prediction/{seq}/{n_fragments}.log",
    conda:
        "envs/pystats.yaml"
    script:
        "scripts/plot_lionelmssq_prediction.py"