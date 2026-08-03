include: "common.smk"


rule plot_prediction:
    input:
        pred_fragments="results/prediction/{modus}/{seq}/{num_replicates}.tsv",
        pred_seq="results/prediction/{modus}/{seq}/{num_replicates}.fasta",
        meta="data/experiment/{seq}/{num_replicates}.preprocessed.meta.yaml",
        sim="data/{modus}/{seq}/{num_replicates}.tsv",
        alphabet="workflow/resources/masses.including_synthetic.tsv",
    output:
        any="results/plots/prediction/{modus}/{seq}/{num_replicates}.html",
        start="results/plots/prediction/{modus}/{seq}/{num_replicates}.start.html",
        end="results/plots/prediction/{modus}/{seq}/{num_replicates}.end.html",
        internal="results/plots/prediction/{modus}/{seq}/{num_replicates}.internal.html",
        all=report(
            "results/plots/prediction/{modus}/{seq}/{num_replicates}.all.html",
            htmlindex="index.html",
            category="Quality control",
            subcategory="{modus} data",
            labels={"sequence": "{seq}", "type": "alignment"},
            caption="../report/quality_control.alignment.rst",
        ),
    log:
        "logs/plots/prediction/{modus}/{seq}/{num_replicates}.log",
    benchmark:
        "benchmarks/plots/prediction/{modus}/{seq}/{num_replicates}.tsv"
    conda:
        "../envs/spectrseqtools.yaml"
    threads: 1
    shell:
        "spectrseqtools plot-fragments "
        "--fragments {input.pred_fragments} "
        "--prediction {input.pred_seq} "
        "--meta {input.meta} "
        "--mixed-fragment-plot {output.any} "
        "--start-fragment-plot {output.start} "
        "--end-fragment-plot {output.end} "
        "--internal-fragment-plot {output.internal} "
        "--combined-plot {output.all} "
        # "--simulation {input.sim} "
        "--alphabet {input.alphabet} "
        "2> {log}"


rule evaluate_custom_simulation:
    input:
        collect_custom_simulations(
            "results/prediction/simulation/{seq}/{num_replicates}.fasta"
        ),
    output:
        "results/evaluation/custom_simulation.tsv",
    log:
        "logs/evaluation/custom_simulation.log",
    benchmark:
        "benchmarks/evaluation/custom_simulation.tsv"
    conda:
        "../envs/spectrseqtools.yaml"
    threads: 1
    script:
        "../scripts/evaluate_prediction.py"


rule plot_evaluation_for_custom_simulation:
    input:
        "results/evaluation/custom_simulation.tsv",
    output:
        donut="results/plots/evaluation/custom_simulation.donut.html",
        bar=report(
            "results/plots/evaluation/custom_simulation.bar.html",
            htmlindex="index.html",
            category="Robustness",
            labels={"type": "simulation data"},
            caption="../report/robustness.data.rst",
        ),
    log:
        "logs/plots/evaluation/custom_simulation.log",
    benchmark:
        "benchmarks/plots/evaluation/custom_simulation.tsv"
    conda:
        "../envs/spectrseqtools.yaml"
    threads: 1
    script:
        "../scripts/plot_evaluation.py"


rule evaluate_comparison_study:
    input:
        lambda wildcards: collect_comparison_studies(
            wildcards.parameter,
            "results/comparison_study/{parameter}/{value}/{seq}/sample.fasta",
        ),
    output:
        "results/comparison_study/{parameter}/evaluation.tsv",
    params:
        mode="simulation",
    log:
        "logs/comparison_study/{parameter}/evaluation.log",
    benchmark:
        "benchmarks/comparison_study/{parameter}/evaluation.log"
    conda:
        "../envs/spectrseqtools.yaml"
    threads: 1
    script:
        "../scripts/evaluate_prediction.py"


rule plot_evaluation_for_comparison_study:
    input:
        "results/comparison_study/{parameter}/evaluation.tsv",
    output:
        donut="results/plots/comparison_study/{parameter}.donut.html",
        bar=report(
            "results/plots/comparison_study/{parameter}.bar.html",
            htmlindex="index.html",
            category="Robustness",
            labels={"type": "comparison", "parameter": "{parameter}"},
            caption="../report/robustness.comparison.rst",
        ),
    log:
        "logs/plots/comparison_study/{parameter}.log",
    benchmark:
        "benchmarks/plots/comparison_study/{parameter}.tsv"
    conda:
        "../envs/spectrseqtools.yaml"
    threads: 1
    script:
        "../scripts/plot_evaluation.py"


rule evaluate_optimization_study:
    input:
        lambda wildcards: collect_optimizations(
            wildcards.parameter,
            "results/optimization/{parameter}/{value}/{seq}/sample.fasta",
        ),
    output:
        "results/optimization/{parameter}/evaluation.tsv",
    params:
        mode="optimization",
    log:
        "logs/optimization/{parameter}/evaluation.log",
    benchmark:
        "benchmarks/optimization/{parameter}/evaluation.log"
    conda:
        "../envs/spectrseqtools.yaml"
    threads: 1
    script:
        "../scripts/evaluate_prediction.py"


rule plot_evaluation_for_optimization_study:
    input:
        "results/optimization/{parameter}/evaluation.tsv",
    output:
        donut="results/plots/optimization/{parameter}.donut.html",
        bar=report(
            "results/plots/optimization/{parameter}.bar.html",
            htmlindex="index.html",
            category="Robustness",
            labels={"type": "optimization", "parameter": "{parameter}"},
            caption="../report/robustness.optimization.rst",
        ),
    log:
        "logs/plots/optimization/{parameter}.log",
    benchmark:
        "benchmarks/plots/optimization/{parameter}.tsv"
    conda:
        "../envs/spectrseqtools.yaml"
    threads: 1
    script:
        "../scripts/plot_evaluation.py"


rule evaluate_random_simulation:
    input:
        collect_random_simulations(
            "results/prediction/simulation/{seq}/{num_replicates}.fasta"
        ),
    output:
        "results/evaluation/random_simulation.tsv",
    log:
        "logs/evaluation/random_simulation.log",
    benchmark:
        "benchmarks/evaluation/random_simulation.tsv"
    conda:
        "../envs/spectrseqtools.yaml"
    threads: 1
    script:
        "../scripts/evaluate_prediction.py"


rule plot_evaluation_for_random_simulation:
    input:
        "results/evaluation/random_simulation.tsv",
    output:
        donut="results/plots/evaluation/random_simulation.donut.html",
        bar="results/plots/evaluation/random_simulation.bar.html",
    log:
        "logs/plots/evaluation/random_simulation.log",
    benchmark:
        "benchmarks/plots/evaluation/random_simulation.tsv"
    conda:
        "../envs/spectrseqtools.yaml"
    threads: 1
    script:
        "../scripts/plot_evaluation.py"


rule evaluate_experiment:
    input:
        collect_experiments(
            "results/prediction/experiment/{seq}/{num_replicates}.fasta"
        ),
    output:
        "results/evaluation/experiment.tsv",
    log:
        "logs/evaluation/experiment.log",
    benchmark:
        "benchmarks/evaluation/experiment.tsv"
    conda:
        "../envs/spectrseqtools.yaml"
    threads: 1
    script:
        "../scripts/evaluate_prediction.py"


rule plot_evaluation_for_experiment:
    input:
        "results/evaluation/experiment.tsv",
    output:
        donut="results/plots/evaluation/experiment.donut.html",
        bar=report(
            "results/plots/evaluation/experiment.bar.html",
            htmlindex="index.html",
            category="Robustness",
            labels={"type": "experiment data"},
            caption="../report/robustness.data.rst",
        ),
    log:
        "logs/plots/evaluation/experiment.log",
    benchmark:
        "benchmarks/plots/evaluation/experiment.tsv"
    conda:
        "../envs/spectrseqtools.yaml"
    threads: 1
    script:
        "../scripts/plot_evaluation.py"


rule plot_spectra:
    input:
        raw_fragments="data/{modus}/{seq}/{num_replicates}.standard_unit_fragments.tsv",
        pred_fragments="results/prediction/{modus}/{seq}/{num_replicates}.tsv",
    output:
        report(
            "results/plots/spectra/{modus}/{seq}/{num_replicates}.html",
            htmlindex="index.html",
            category="Quality control",
            subcategory="{modus} data",
            labels={
                "sequence": "{seq}",
                "type": "spectrum",
            },
            caption="../report/quality_control.spectrum.rst",
        ),
    log:
        "logs/plots/spectra/{modus}/{seq}/{num_replicates}.log",
    benchmark:
        "benchmarks/plots/spectra/{modus}/{seq}/{num_replicates}.tsv"
    conda:
        "../envs/spectrseqtools.yaml"
    threads: 1
    script:
        "../scripts/plot_spectra.py"


rule plot_singletons:
    input:
        raw_data="data/experiment/{seq}/{num_replicates}.raw",
        meta="data/experiment/{seq}/{num_replicates}.meta.yaml",
        alphabet="workflow/resources/masses.including_synthetic.tsv",
        singletons="data/experiment/{seq}/{num_replicates}.singletons.tsv",
    output:
        all=report(
            "results/plots/singletons/{seq}/{num_replicates}.html",
            htmlindex="index.html",
            category="Quality control",
            subcategory="experiment data",
            labels={
                "sequence": "{seq}",
                "type": "singletons",
            },
            caption="../report/quality_control.singletons.rst",
        ),
        single="results/plots/singletons/{seq}/{num_replicates}_single/scan_0.html",
    params:
        scan_dir=subpath(output.single, parent=True),
    log:
        "logs/plots/singletons/{seq}/{num_replicates}.log",
    benchmark:
        "benchmarks/plots/singletons/{seq}/{num_replicates}.tsv"
    conda:
        "../envs/spectrseqtools.yaml"
    threads: 1
    shell:
        "spectrseqtools plot-singletons "
        "--input {input.raw_data} "
        "--meta {input.meta} "
        "--scan-dir {params.scan_dir} "
        "--output-path {output.all} "
        "--alphabet {input.alphabet} "
        "2> {log}"


rule evaluate_run_statistics_for_simulations:
    input:
        benchmarks=collect_comparison_studies(
            "num_replicates",
            "benchmarks/comparison_study/num_replicates/{value}/{seq}/prediction.tsv",
        ),
        fragments=collect_comparison_studies(
            "num_replicates",
            "comparison_study/num_replicates/{value}/{seq}/sample.tsv",
        ),
    output:
        "results/comparison_study/stats.tsv",
    log:
        "logs/comparison_study/stats.log",
    benchmark:
        "benchmarks/comparison_study/stats.log"
    conda:
        "../envs/spectrseqtools.yaml"
    threads: 1
    script:
        "../scripts/evaluate_run_statistics.py"


rule evaluate_run_statistics_for_experiments:
    input:
        benchmarks=collect_experiments(
            "benchmarks/prediction/experiment/{seq}/{num_replicates}.tsv"
        ),
        fragments=collect_experiments(
            "data/experiment/{seq}/{num_replicates}.tsv",
        ),
    output:
        "results/comparison_study/experiments.stats.tsv",
    log:
        "logs/comparison_study/experiments.stats.log",
    benchmark:
        "benchmarks/comparison_study/experiments.stats.log"
    conda:
        "../envs/spectrseqtools.yaml"
    threads: 1
    script:
        "../scripts/evaluate_run_statistics.py"


rule plot_runtime:
    input:
        sim="results/comparison_study/stats.tsv",
        exp="results/comparison_study/experiments.stats.tsv",
    output:
        report(
            "results/plots/evaluation/runtime.html",
            htmlindex="index.html",
            category="Robustness",
            labels={"type": "runtime"},
            # caption="../report/robustness.data.rst",
        ),
    params:
        mode="runtime",
    log:
        "logs/plots/evaluation/runtime.log",
    benchmark:
        "benchmarks/plots/evaluation/runtime.tsv"
    conda:
        "../envs/spectrseqtools.yaml"
    threads: 1
    script:
        "../scripts/plot_run_statistics.py"


rule plot_memory:
    input:
        sim="results/comparison_study/stats.tsv",
        exp="results/comparison_study/experiments.stats.tsv",
    output:
        report(
            "results/plots/evaluation/memory.html",
            htmlindex="index.html",
            category="Robustness",
            labels={"type": "memory"},
            # caption="../report/robustness.data.rst",
        ),
    params:
        mode="memory",
    log:
        "logs/plots/evaluation/memory.log",
    benchmark:
        "benchmarks/plots/evaluation/memory.tsv"
    conda:
        "../envs/spectrseqtools.yaml"
    threads: 1
    script:
        "../scripts/plot_run_statistics.py"
