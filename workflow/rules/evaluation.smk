rule plot_lionelmssq_prediction:
    input:
        pred_fragments="results/lionelmssq/{seq}/{n_fragments}.tsv",
        pred_seq="results/lionelmssq/{seq}/{n_fragments}.fasta",
        sim="data/simulation/{seq}/{n_fragments}.tsv",
    output:
        "results/plots/lionelmssq_prediction/{seq}/{n_fragments}.html",
    log:
        "logs/plots/lionelmssq_prediction/{seq}/{n_fragments}.log",
    conda:
        "../envs/lionelmssq.yaml"
    script:
        "../scripts/plot_lionelmssq_prediction.py"
