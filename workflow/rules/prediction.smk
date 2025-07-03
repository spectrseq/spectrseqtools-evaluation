rule lionelmssq:
    input:
        fragments="data/{modus}/{seq}/{n_fragments}.tsv",
    output:
        fragment_predictions="results/lionelmssq/{modus}/{seq}/{n_fragments}.tsv",
        sequence="results/lionelmssq/{modus}/{seq}/{n_fragments}.fasta",
    params:
        seq_len=get_seq_len,
    log:
        "logs/lionelmssq/{modus}/{seq}/{n_fragments}.log",
    conda:
        "../envs/lionelmssq.yaml"
    threads: 64
    shell:
        "lionelmssq --fragments {input.fragments} --seq-len {params.seq_len} "
        "--fragment-predictions {output.fragment_predictions} "
        "--sequence-prediction {output.sequence} "
        "--sequence-name 'lionelmssq_prediction_from_sim_{wildcards.seq}' "
        "--solver cbc "
        "--threads {threads} "
        "2> {log}"
