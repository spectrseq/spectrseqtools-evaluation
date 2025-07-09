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


rule plot_simulated_fragments:
    input:
        config=workflow.source_path("../resources/datavzrd/simulation.yaml"),
        simulation="data/simulation/{seq}/{n_fragments}.tsv",
    output:
        report(
            directory("results/plots/simulated_fragments/{seq}/{n_fragments}"),
            htmlindex="index.html",
            category="Simulation",
            labels={"seq": "{seq}", "n_fragments": "{n_fragments}"},
        ),
    log:
        "logs/plots/simulated_fragments/{seq}/{n_fragments}.log",
    wrapper:
        "v7.2.0/utils/datavzrd"
