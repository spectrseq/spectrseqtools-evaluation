rule simulate_measurement:
    input:
        workflow.source_path("../resources/masses_all.tsv"),
    output:
        "results/simulation/{seq}/{n_fragments}.tsv",
    log:
        "logs/simulation/{seq}/{n_fragments}.log",
    conda:
        "../envs/pystats.yaml"
    script:
        "../scripts/simulate_rna_measurement.py"


rule datavzrd:
    input:
        config=workflow.source_path("../resources/datavzrd/simulation.yaml"),
        simulation="results/simulation/{seq}/{n_fragments}.tsv",
    output:
        report(
            directory("results/datavzrd/simulation/{seq}/{n_fragments}"),
            htmlindex="index.html",
            category="Simulation",
            labels={"seq": "{seq}", "n_fragments": "{n_fragments}"},
        ),
    log:
        "logs/datavzrd/simulation/{seq}/{n_fragments}.log",
    wrapper:
        "v5.0.2/utils/datavzrd"
        #"v5.0.1-1-g34a454a/utils/datavzrd"
