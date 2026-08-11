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


def generate_random_sequence_and_seed_pair(
    seq_len, modification_rate=0, modifications=None
):
    # If no modifications are given, simulate over unmodified bases only
    if modifications is None:
        return "".join(
            [random.choice(["A", "U", "G", "C"]) for _ in range(seq_len)]
        ), random.choice(range(10000))

    # Read modifications from file
    with open(modifications, "r") as file:
        lines = file.readlines()[1:]
    modified_nucleoside_names = [
        line.split("\t")[0]
        for line in lines
        if line.split("\t")[0] not in ["U", "A", "G", "C"]
    ]

    # Define probabilities for different bases (using modification rate)
    weights = [(1.0 - modification_rate) / 4] * 4 + [
        modification_rate / len(modified_nucleoside_names)
    ] * len(modified_nucleoside_names)

    return "".join(
        random.choices(
            ["A", "U", "G", "C"] + modified_nucleoside_names,
            weights=weights,
            k=seq_len,
        )
    ), random.choice(range(10000))


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
    retval = []

    if lookup(dpath="comparison/studies", within=config) is None:
        print("No random simulation data given.")
        return retval

    if param not in lookup(dpath="comparison/studies", within=config):
        print(f"No data for {param} given.")
        return retval

    if lookup(dpath=f"comparison/studies/{param}", within=config) is None:
        print("No comparison studies given.")
        return []

    sim = lookup(dpath=f"comparison/studies/{param}", within=config)
    for modification_rate in sim["modification_rate"]:
        for seq_len in sim["sequence_length"]:
            random.seed(lookup(dpath="comparison/seed", within=config))
            sequences = [
                generate_random_sequence_and_seed_pair(
                    seq_len=random.choice(range(10, 21)) if seq_len == -1 else seq_len,
                    modification_rate=modification_rate,
                    modifications=workflow.source_path("../resources/masses.tsv"),
                )
                for _ in range(lookup(dpath="comparison/num_sequences", within=config))
            ]

            for seq in sequences:
                if param == "modification_rate":
                    values = [modification_rate]
                elif param == "sequence_length":
                    values = [seq_len]
                else:
                    values = sim[param]
                for value in values:
                    file_name = f"comparison_study/{param}/{value}/{seq[0]}/seed.txt"
                    os.makedirs(os.path.dirname(file_name), exist_ok=True)
                    with open(file_name, "w") as f:
                        f.write(str(seq[1]))

                    retval += collect(
                        patterns,
                        parameter=param,
                        value=value,
                        seq=seq[0],
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
