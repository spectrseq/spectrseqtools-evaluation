import random
import os


def has_custom_percentile(wildcards):
    for item in lookup(dpath="experiment", within=config):
        if item["seq"] == wildcards.seq:
            return "percentile" in item.keys()
    return False


def get_custom_percentile(wildcards):
    for item in lookup(dpath="experiment", within=config):
        if item["seq"] == wildcards.seq:
            return item["percentile"]
    return None


def collect_simulations(*patterns):
    if lookup(dpath="simulation", within=config) is None:
        print("No custom simulation data given.")
        return []

    retval = []
    for item in lookup(dpath="simulation", within=config):
        retval += collect(
            patterns,
            # seq=lookup(dpath=f"simulation/{item}/seq", within=config),
            seq=item,
            num_replicates="sample",
        )
    return retval


def collect_comparison_studies(param: str, *patterns):
    if lookup(dpath="comparison/studies", within=config) is None:
        print("No random simulation data given.")
        return []

    if param not in lookup(dpath="comparison/studies", within=config):
        print(f"No data for {param} given.")
        return []

    if lookup(dpath=f"comparison/studies/{param}", within=config) is None:
        print("No comparison studies given.")
        return []

    retval = []
    sim = lookup(dpath=f"comparison/studies/{param}", within=config)
    for modification_rate in sim["modification_rate"]:
        for seq_len in sim["sequence_length"]:
            for id in range(lookup(dpath="comparison/num_sequences", within=config)):
                if param == "modification_rate":
                    values = [modification_rate]
                elif param == "sequence_length":
                    values = [seq_len]
                else:
                    values = sim[param]

                for value in values:
                    retval += collect(
                        patterns,
                        parameter=param,
                        value=value,
                        id=id + 1,
                    )
    return retval


def collect_experiments(*patterns):
    if lookup(dpath="experiment", within=config) is None:
        print("No experimental data given.")
        return []

    return [
        collect(
            patterns,
            seq=item["seq"],
            num_replicates=item["fragments"],
        )
        for item in lookup(dpath="experiment", within=config)
    ]


def collect_optimizations(param: str, *patterns):
    if lookup(dpath="experiment", within=config) is None:
        print("No evaluation data given.")
        return []

    if param not in lookup(dpath="optimization", within=config):
        print(f"No data for {param} given.")
        return []

    if lookup(dpath=f"optimization/{param}", within=config) is None:
        print("No optimization studies given.")
        return []

    retval = []
    for value in lookup(dpath=f"optimization/{param}/{param}", within=config):
        for item in lookup(dpath="experiment", within=config):
            retval += collect(
                patterns,
                parameter=param,
                value=value,
                seq=item["seq"],
                num_replicates=item["fragments"],
            )

    return retval
