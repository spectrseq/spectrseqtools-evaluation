rule lionelmssq:
    input:
        "results/simulation/{seq}/{n_fragments}.tsv",
    output:
        fragment_predictions="results/lionelmssq/{seq}/{n_fragments}.tsv",
        sequence="results/lionelmssq/{seq}/{n_fragments}.fasta",
    params:
        seq_len=get_seq_len,
    log:
        "logs/lionelmssq/{seq}/{n_fragments}.log",
    threads: 64
    shell:
        "lionelmssq --fragments {input} --seq-len {params.seq_len} "
        "--fragment-predictions {output.fragment_predictions} "
        "--sequence-prediction {output.sequence} "
        "--sequence-name 'lionelmssq_prediction_from_sim_{wildcards.seq}' "
        "--solver gurobi "
        "--threads {threads} "
        "2> {log}"
