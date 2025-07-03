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


rule plot_fragment_data:
    input:
        config=workflow.source_path("../resources/datavzrd/simulation.yaml"),
        simulation="data/{modus}/{seq}/{n_fragments}.tsv",
    output:
        report(
            directory("results/plots/fragments/{modus}/{seq}/{n_fragments}"),
            htmlindex="index.html",
            category="{modus}",
            labels={"seq": "{seq}", "n_fragments": "{n_fragments}"},
        ),
    log:
        "logs/plots/fragments/{modus}/{seq}/{n_fragments}.log",
    wrapper:
        "v5.0.2/utils/datavzrd"
