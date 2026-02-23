rule simulate_for_comparison_study:
    input:
        nucleosides=workflow.source_path("../resources/masses.tsv"),
        elements=workflow.source_path("../resources/element_masses.tsv"),
    output:
        fragments="comparison_study/{parameter}/{value}/{seq}/sample.tsv",
        singletons="comparison_study/{parameter}/{value}/{seq}/sample.singletons.tsv",
        meta="comparison_study/{parameter}/{value}/{seq}/sample.meta.yaml",
    params:
        dir="comparison_study/{parameter}/{value}/{seq}",
        num_copies=lambda wildcards: wildcards.value
        if wildcards.parameter == "num_copies"
        else lookup(
            dpath=f"comparison/studies/{wildcards.parameter}/num_copies",
            within=config,
        )[0],
        max_singletons=lambda wildcards: wildcards.value
        if wildcards.parameter == "max_singletons"
        else lookup(
            dpath=f"comparison/studies/{wildcards.parameter}/max_singletons",
            within=config,
        )[0],
        ghost_rate=lambda wildcards: wildcards.value
        if wildcards.parameter == "ghost_rate"
        else lookup(
            dpath=f"comparison/studies/{wildcards.parameter}/ghost_rate",
            within=config,
        )[0],
        rel_error_rate=lambda wildcards: wildcards.value
        if wildcards.parameter == "rel_error_rate"
        else lookup(
            dpath=f"comparison/studies/{wildcards.parameter}/rel_error_rate",
            within=config,
        )[0],
    log:
        "logs/comparison_study/{parameter}/{value}/{seq}/simulation.log",
    benchmark:
        "benchmarks/comparison_study/{parameter}/{value}/{seq}/simulation.tsv"
    conda:
        "../envs/spectrseqtools.yaml"
    threads: 1
    script:
        "../scripts/simulate_fragments.py"


rule simulate_measurement:
    input:
        nucleosides=workflow.source_path("../resources/masses.tsv"),
        elements=workflow.source_path("../resources/element_masses.tsv"),
    output:
        fragments="data/simulation/{seq}/{n_fragments}.tsv",
        singletons="data/simulation/{seq}/{n_fragments}.singletons.tsv",
        meta="data/simulation/{seq}/{n_fragments}.meta.yaml",
    params:
        dir=None,
        num_copies=lambda wildcards: wildcards.n_fragments,
        max_singletons=lookup(
            dpath="fragmentation_params/max_singletons",
            within=config,
        ),
        ghost_rate=lookup(
            dpath="fragmentation_params/ghost_rate",
            within=config,
        ),
        rel_error_rate=lookup(
            dpath="fragmentation_params/rel_error_rate",
            within=config,
        ),
    log:
        "logs/simulation/{seq}/{n_fragments}.log",
    benchmark:
        "benchmarks/simulation/{seq}/{n_fragments}.tsv"
    conda:
        "../envs/spectrseqtools.yaml"
    threads: 1
    script:
        "../scripts/simulate_fragments.py"


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
        "benchmarks/plots/simulated_fragments/{seq}/{n_fragments}.tsv"
    threads: 1
    wrapper:
        "v7.2.0/utils/datavzrd"
