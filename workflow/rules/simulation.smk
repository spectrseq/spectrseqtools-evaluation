rule simulate_measurement:
    input:
        nucleosides=workflow.source_path("../resources/masses.tsv"),
        elements=workflow.source_path("../resources/element_masses.tsv"),
    output:
        fragments="data/simulation/{seq}/{n_fragments}.tsv",
        meta="data/simulation/{seq}/{n_fragments}.meta.yaml",
    log:
        "logs/simulation/{seq}/{n_fragments}.log",
    benchmark:
        "benchmarks/simulation/{seq}/{n_fragments}.tsv",
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
    benchmark:
        "benchmarks/plots/simulated_fragments/{seq}/{n_fragments}.tsv",
    wrapper:
        "v7.2.0/utils/datavzrd"
