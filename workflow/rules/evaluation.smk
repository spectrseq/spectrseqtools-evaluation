rule plot_lionelmssq_prediction:
    input:
        pred_fragments="results/lionelmssq/{modus}/{seq}/{n_fragments}.tsv",
        pred_seq="results/lionelmssq/{modus}/{seq}/{n_fragments}.fasta",
        sim="data/{modus}/{seq}/{n_fragments}.tsv",
    output:
        "results/plots/lionelmssq_prediction/{modus}/{seq}/{n_fragments}.html",
    log:
        "logs/plots/lionelmssq_prediction/{modus}/{seq}/{n_fragments}.log",
    conda:
        "../envs/lionelmssq.yaml"
    script:
        "../scripts/plot_lionelmssq_prediction.py"
