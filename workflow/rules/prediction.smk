rule prediction_experiment:
    input:
        fragments="data/experiment/{seq}/{n_fragments}.raw",
        meta="data/experiment/{seq}/{n_fragments}.meta.yaml",
    output:
        fragments="data/experiment/{seq}/{n_fragments}.tsv",
        singletons="data/experiment/{seq}/{n_fragments}.singletons.tsv",
        su_fragments="data/experiment/{seq}/{n_fragments}.standard_unit_fragments.tsv",
        predictions="results/prediction/experiment/{seq}/{n_fragments}.tsv",
        sequence="results/prediction/experiment/{seq}/{n_fragments}.fasta",
    params:
        solver=config["solver"],
    log:
        "logs/prediction/experiment/{seq}/{n_fragments}.log",
    benchmark:
        "benchmarks/prediction/experiment/{seq}/{n_fragments}.tsv"
    conda:
        "../envs/spectrseqtools.yaml"
    threads: 1
    shell:
        "spectrseqtools --fragments {input.fragments} --meta {input.meta} "
        "--fragment-predictions {output.predictions} "
        "--sequence-prediction {output.sequence} "
        "--sequence-name 'spectrseqtools_prediction_{wildcards.seq}' "
        "--solver {params.solver} "
        "--threads {threads} "
        "2> {log}"


rule prediction_simulation:
    input:
        fragments="data/simulation/{seq}/{n_fragments}.tsv",
        meta="data/simulation/{seq}/{n_fragments}.meta.yaml",
        singletons="data/simulation/{seq}/{n_fragments}.singletons.tsv",
    output:
        su_fragments="data/simulation/{seq}/{n_fragments}.standard_unit_fragments.tsv",
        predictions="results/prediction/simulation/{seq}/{n_fragments}.tsv",
        sequence="results/prediction/simulation/{seq}/{n_fragments}.fasta",
    params:
        solver=config["solver"],
    log:
        "logs/prediction/simulation/{seq}/{n_fragments}.log",
    benchmark:
        "benchmarks/prediction/simulation/{seq}/{n_fragments}.tsv"
    conda:
        "../envs/spectrseqtools.yaml"
    threads: 1
    shell:
        "spectrseqtools --fragments {input.fragments} --meta {input.meta} "
        "--fragment-predictions {output.predictions} "
        "--sequence-prediction {output.sequence} "
        "--sequence-name 'spectrseqtools_prediction_from_sim_{wildcards.seq}' "
        "--solver {params.solver} "
        "--threads {threads} "
        "2> {log}"


rule prediction_comparison_study:
    input:
        fragments="comparison_study/{parameter}/{value}/{seq}/sample.tsv",
        meta="comparison_study/{parameter}/{value}/{seq}/sample.meta.yaml",
    output:
        predictions="results/comparison_study/{parameter}/{value}/{seq}/sample.tsv",
        sequence="results/comparison_study/{parameter}/{value}/{seq}/sample.fasta",
    params:
        solver=config["solver"],
    log:
        "logs/comparison_study/{parameter}/{value}/{seq}/prediction.log",
    benchmark:
        "benchmarks/comparison_study/{parameter}/{value}/{seq}/prediction.tsv"
    conda:
        "../envs/spectrseqtools.yaml"
    threads: 1
    shell:
        "spectrseqtools --fragments {input.fragments} --meta {input.meta} "
        "--fragment-predictions {output.predictions} "
        "--sequence-prediction {output.sequence} "
        "--sequence-name 'spectrseqtools_prediction_from_sim_{wildcards.seq}' "
        "--solver {params.solver} "
        "--threads {threads} "
        "2> {log}"


rule prediction_optimization_study:
    input:
        fragments="data/experiment/{seq}/0.raw",
        meta="data/experiment/{seq}/0.meta.yaml",
    output:
        meta="results/optimization/{parameter}/{value}/{seq}/0.preprocessed.meta.yaml",
        fragments="results/optimization/{parameter}/{value}/{seq}/0.tsv",
        singletons="results/optimization/{parameter}/{value}/{seq}/0.singletons.tsv",
        predictions="results/optimization/{parameter}/{value}/{seq}/sample.tsv",
        sequence="results/optimization/{parameter}/{value}/{seq}/sample.fasta",
    params:
        solver=config["solver"],
        intensity_cutoff=lambda wildcards: wildcards.value
        if wildcards.parameter == "intensity_cutoff"
        else lookup(
            dpath=f"optimization/{wildcards.parameter}/intensity_cutoff",
            within=config,
        )[0],
        lp_timeout_long=lambda wildcards: wildcards.value
        if wildcards.parameter == "lp_timeout_long"
        else lookup(
            dpath=f"optimization/{wildcards.parameter}/lp_timeout_long",
            within=config,
        )[0],
        lp_timeout_short=lambda wildcards: wildcards.value
        if wildcards.parameter == "lp_timeout_short"
        else lookup(
            dpath=f"optimization/{wildcards.parameter}/lp_timeout_short",
            within=config,
        )[0],
        dir=subpath(output.fragments, parent=True),
        # dir="results/optimization/{parameter}/{value}/{seq}/",
    log:
        "logs/optimization/{parameter}/{value}/{seq}/sample.log",
    benchmark:
        "benchmarks/optimization/{parameter}/{value}/{seq}/sample.tsv"
    conda:
        "../envs/spectrseqtools.yaml"
    threads: 1
    shell:
        "spectrseqtools --fragments {input.fragments} --meta {input.meta} "
        "--fragment-predictions {output.predictions} "
        "--sequence-prediction {output.sequence} "
        "--sequence-name 'spectrseqtools_prediction_{wildcards.seq}' "
        "--solver {params.solver} "
        "--threads {threads} "
        "--cutoff-percentile {params.intensity_cutoff} "
        "--lp-timeout-long {params.lp_timeout_long} "
        "--lp-timeout-short {params.lp_timeout_short} "
        "--output-dir {params.dir} "
        "2> {log}"
