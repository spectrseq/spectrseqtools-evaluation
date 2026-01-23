import random
import os


wildcard_constraints:
    n_fragments="[0-9]+",


def generate_random_sequence_and_seed_pair(
    seq_len, mutation_rate=0, modifications=None
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

    # Define probabilities for different bases (using mutation rate)
    weights = [(1.0 - mutation_rate) / 4] * 4 + [
        mutation_rate / len(modified_nucleoside_names)
    ] * len(modified_nucleoside_names)

    return "".join(
        random.choices(
            ["A", "U", "G", "C"] + modified_nucleoside_names,
            weights=weights,
            k=seq_len,
        )
    ), random.choice(range(10000))


def generate_random_sequences(seq_len, mutation_rate=0, modified_nucleosides=None):
    if modified_nucleosides is not None:
        with open(modified_nucleosides, "r") as file:
            lines = file.readlines()[1:]
        modified_nucleoside_names = [
            line.split("\t")[0]
            for line in lines
            if line.split("\t")[0] not in ["U", "A", "G", "C"]
        ]

        # Define probabilities for different bases (using mutation rate)
        weights = [(1.0 - mutation_rate) / 4] * 4 + [
            mutation_rate / len(modified_nucleoside_names)
        ] * len(modified_nucleoside_names)

        return "".join(
            random.choices(
                ["A", "U", "G", "C"] + modified_nucleoside_names,
                weights=weights,
                k=seq_len,
            )
        )

    else:
        return "".join([random.choice(["A", "U", "G", "C"]) for _ in range(seq_len)])


def collect_custom_simulations(*patterns):
    if lookup(dpath="simulation/custom", within=config) is None:
        print("No custom simulation data given.")
        return []

    return [
        collect(
            patterns,
            seq=item["seq"],
            n_fragments=n_fragments,
        )
        for item in lookup(dpath="simulation/custom", within=config)
        for n_fragments in item["n_fragments"]
    ]


def collect_random_simulations(*patterns):
    retval = []

    if lookup(dpath="simulation/random/strata", within=config) is None:
        print("No random simulation data given.")
        return retval

    random.seed(lookup(dpath="simulation/random", within=config)["seed"])
    for sim in lookup(dpath="simulation/random/strata", within=config):
        sequences = [
            generate_random_sequence_and_seed_pair(
                seq_len=random.choice(range(10, 21)),
                mutation_rate=lookup(dpath="simulation/random", within=config)[
                    "mutation_rate"
                ],
                modifications=workflow.source_path("../resources/masses.tsv"),
            )
            for _ in range(sim["num_seq"])
        ]

        for seq in sequences:
            file_name = f"data/simulation/{seq[0]}/seed.txt"
            os.makedirs(os.path.dirname(file_name), exist_ok=True)
            with open(file_name, "w") as f:
                f.write(str(seq[1]))

            retval += [
                collect(
                    patterns,
                    seq=seq[0],
                    n_fragments=n_fragments,
                )
                for n_fragments in sim["n_fragments"]
            ]

    return retval


def collect_comparison_studies(param: str, *patterns):
    retval = []

    if lookup(dpath="comparison/studies", within=config) is None:
        print("No random simulation data given.")
        return retval

    sim = lookup(dpath=f"comparison/studies/{param}", within=config)
    for mutation_rate in sim["mutation_rate"]:
        random.seed(lookup(dpath="comparison/seed", within=config))
        sequences = [
            generate_random_sequence_and_seed_pair(
                seq_len=random.choice(range(10, 21)),
                mutation_rate=mutation_rate,
                modifications=workflow.source_path("../resources/masses.tsv"),
            )
            for _ in range(lookup(dpath="comparison/num_sequences", within=config))
        ]

        for seq in sequences:
            values = sim[param] if param != "mutation_rate" else [mutation_rate]
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


def collect_simulations(*patterns):
    retval = []

    retval += collect_custom_simulations(*patterns)
    retval += collect_random_simulations(*patterns)

    return retval


def collect_experiments(*patterns):
    if lookup(dpath="experiment", within=config) is None:
        print("No experimental data given.")
        return []

    return [
        collect(
            patterns,
            seq=item["seq"],
            n_fragments=item["fragments"],
        )
        for item in lookup(dpath="experiment", within=config)
    ]
