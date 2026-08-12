wildcard_constraints:
    value="[0-9]+",


rule simulate_metadata_for_comparison_study:
    input:
        alphabet=workflow.source_path("../resources/masses.tsv"),
    output:
        dir=directory("comparison_study/{parameter}/{value}/"),
        meta=expand("comparison_study/{{parameter}}/{{value}}/sim_{id}/sample.meta.yaml",
        id=range(1, lookup(
            dpath="comparison/num_sequences",
            within=config,
        )+1)),
    log:
        "logs/comparison_study/{parameter}/{value}/metadata_simulation.log",
    benchmark:
        "benchmarks/comparison_study/{parameter}/{value}/metadata_simulation.tsv"
    conda:
        "../envs/spectrseqtools.yaml"
    threads: 1
    params:
        start_tag=config["fragmentation_params"]["5_prime_tag"],
        end_tag=config["fragmentation_params"]["3_prime_tag"],
        num_seqs=lookup(
            dpath="comparison/num_sequences",
            within=config,
        ),
        seq_len=lambda wildcards: (
            wildcards.value
            if wildcards.parameter == "sequence_length"
            else lookup(
                dpath=f"comparison/studies/{wildcards.parameter}/sequence_length",
                within=config,
            )[0]
        ),
        mod_rate=lambda wildcards: (
            wildcards.value
            if wildcards.parameter == "modification_rate"
            else lookup(
                dpath=f"comparison/studies/{wildcards.parameter}/modification_rate",
                within=config,
            )[0]
        ),
        seed=lookup(
            dpath="comparison/seed",
            within=config,
        ),
    shell:
        "spectrseqtools simulation random "
        "--num-sequences {params.num_seqs} "
        "--output-dir {output.dir} "
        "--start-tag {params.start_tag} "
        "--end-tag {params.end_tag} "
        "--sequence-length {params.seq_len} "
        "--modification-rate {params.mod_rate} "
        "--alphabet {input.alphabet} "
        "--global-seed {params.seed} "
        "2> {log}"


rule simulate_for_comparison_study:
    input:
        elements=workflow.source_path("../resources/element_masses.tsv"),
        meta="comparison_study/{parameter}/{value}/{seq}/sample.meta.yaml",
    output:
        fragments="comparison_study/{parameter}/{value}/{seq}/sample.tsv",
        singletons="comparison_study/{parameter}/{value}/{seq}/sample.singletons.tsv",
        meta="comparison_study/{parameter}/{value}/{seq}/sample.preprocessed.meta.yaml",
    log:
        "logs/comparison_study/{parameter}/{value}/{seq}/simulation.log",
    benchmark:
        "benchmarks/comparison_study/{parameter}/{value}/{seq}/simulation.tsv"
    conda:
        "../envs/spectrseqtools.yaml"
    threads: 1
    params:
        num_replicates=lambda wildcards: (
            wildcards.value
            if wildcards.parameter == "num_replicates"
            else lookup(
                dpath=f"comparison/studies/{wildcards.parameter}/num_replicates",
                within=config,
            )[0]
        ),
        max_singletons=lambda wildcards: (
            wildcards.value
            if wildcards.parameter == "max_singletons"
            else lookup(
                dpath=f"comparison/studies/{wildcards.parameter}/max_singletons",
                within=config,
            )[0]
        ),
        phantom_rate=lambda wildcards: (
            wildcards.value
            if wildcards.parameter == "phantom_rate"
            else lookup(
                dpath=f"comparison/studies/{wildcards.parameter}/phantom_rate",
                within=config,
            )[0]
        ),
        noise_rate=lambda wildcards: (
            wildcards.value
            if wildcards.parameter == "noise_rate"
            else lookup(
                dpath=f"comparison/studies/{wildcards.parameter}/noise_rate",
                within=config,
            )[0]
        ),
        config=lookup(
            dpath="fragmentation_params",
            within=config,
        ),
    shell:
        "spectrseqtools simulation fragments "
        "--elements {input.elements} "
        "--input {input.meta} "
        "--fragments {output.fragments} "
        "--singletons {output.singletons} "
        "--meta {output.meta} "
        "--num-replicates {params.num_replicates} "
        "--max-singletons {params.max_singletons} "
        "--phantom-rate {params.phantom_rate} "
        "--noise-rate {params.noise_rate} "
        '--config "{params.config}" '
        "2> {log}"


rule simulate_metadata_for_simulation:
    input:
        alphabet=workflow.source_path("../resources/masses.tsv"),
    output:
        meta="data/simulation/{seq}/sample.meta.yaml",
    log:
        "logs/simulation/{seq}/metadata_simulation.log",
    benchmark:
        "benchmarks/simulation/{seq}/metadata_simulation.tsv"
    conda:
        "../envs/spectrseqtools.yaml"
    threads: 1
    params:
        dir=subpath(output.meta, parent=True),
        start_tag=config["fragmentation_params"]["5_prime_tag"],
        end_tag=config["fragmentation_params"]["3_prime_tag"],
        seq=lookup(dpath="simulation/{seq}/seq", within=config),
    shell:
        "spectrseqtools simulation custom "
        "--sequence {params.seq} "
        "--output-dir {params.dir} "
        "--start-tag {params.start_tag} "
        "--end-tag {params.end_tag} "
        "2> {log}"


rule simulate_custom_fragments:
    input:
        elements=workflow.source_path("../resources/element_masses.tsv"),
        meta="data/simulation/{seq}/{num_replicates}.meta.yaml",
    output:
        fragments="data/simulation/{seq}/{num_replicates}.tsv",
        singletons="data/simulation/{seq}/{num_replicates}.singletons.tsv",
        meta="data/simulation/{seq}/{num_replicates}.preprocessed.meta.yaml",
    log:
        "logs/simulation/{seq}/{num_replicates}.log",
    benchmark:
        "benchmarks/simulation/{seq}/{num_replicates}.tsv"
    conda:
        "../envs/spectrseqtools.yaml"
    threads: 1
    params:
        num_replicates=lookup(
            dpath="simulation/{seq}/num_replicates",
            within=config,
        ),
        max_singletons=lookup(
            dpath="fragmentation_params/max_singletons",
            within=config,
        ),
        phantom_rate=lookup(
            dpath="fragmentation_params/phantom_rate",
            within=config,
        ),
        noise_rate=lookup(
            dpath="fragmentation_params/noise_rate",
            within=config,
        ),
        config=lookup(
            dpath="fragmentation_params",
            within=config,
        ),
    shell:
        "spectrseqtools simulation fragments "
        "--elements {input.elements} "
        "--input {input.meta} "
        "--fragments {output.fragments} "
        "--singletons {output.singletons} "
        "--meta {output.meta} "
        "--num-replicates {params.num_replicates} "
        "--max-singletons {params.max_singletons} "
        "--phantom-rate {params.phantom_rate} "
        "--noise-rate {params.noise_rate} "
        '--config "{params.config}" '
        "2> {log}"


rule plot_simulated_fragments:
    input:
        config=workflow.source_path("../resources/datavzrd/simulation.yaml"),
        simulation="data/simulation/{seq}/{num_replicates}.tsv",
    output:
        report(
            directory("results/plots/simulated_fragments/{seq}/{num_replicates}"),
            htmlindex="index.html",
            category="Simulation",
            labels={"seq": "{seq}", "num_replicates": "{num_replicates}"},
        ),
    log:
        "logs/plots/simulated_fragments/{seq}/{num_replicates}.log",
    benchmark:
        "benchmarks/plots/simulated_fragments/{seq}/{num_replicates}.tsv"
    threads: 1
    wrapper:
        "v7.2.0/utils/datavzrd"
