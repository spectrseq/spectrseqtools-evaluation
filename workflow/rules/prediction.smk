rule lionelmssq:
    input:
        fragments="results/simulation/{seq}/{n_fragments}.tsv",
        bin=lookup(dpath="dev/envs/lionelmssq", within=config) + "/bin/lionelmssq",
    output:
        fragment_predictions="results/lionelmssq/{seq}/{n_fragments}.tsv",
        sequence="results/lionelmssq/{seq}/{n_fragments}.fasta",
    params:
        seq_len=get_seq_len,
    log:
        "logs/lionelmssq/{seq}/{n_fragments}.log",
    threads: 64
    shell:
        "{input.bin} --fragments {input.fragments} --seq-len {params.seq_len} "
        "--fragment-predictions {output.fragment_predictions} "
        "--sequence-prediction {output.sequence} "
        "--sequence-name 'lionelmssq_prediction_from_sim_{wildcards.seq}' "
        "--solver gurobi "
        "--threads {threads} "
        "2> {log}"
