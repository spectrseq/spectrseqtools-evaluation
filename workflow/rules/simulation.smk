rule simulate_measurement:
    input:
        nucleosides=workflow.source_path("../resources/masses_all.tsv"),
        elements=workflow.source_path("../resources/element_masses.tsv"),
    output:
        "data/simulation/{seq}/{n_fragments}.tsv",
    log:
        "logs/simulation/{seq}/{n_fragments}.log",
    conda:
        "../envs/pystats.yaml"
    script:
        "../scripts/simulate_rna_measurement.py"


rule datavzrd:
    input:
        config=workflow.source_path("../resources/datavzrd/simulation.yaml"),
        simulation="data/simulation/{seq}/{n_fragments}.tsv",
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
