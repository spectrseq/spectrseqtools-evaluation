import re, random


wildcard_constraints:
    n_fragments="[0-9]+",


_NUCLEOSIDE_RE = re.compile(r"\d*[ACGU]")


def get_seq_len(wildcards):
    return len(_NUCLEOSIDE_RE.findall(wildcards.seq))


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


def collect_simulations(*patterns):
    retval = []

    if lookup(dpath="simulation/custom", within=config) is not None:
        retval += [
            collect(
                patterns,
                seq=item["seq"],
                n_fragments=n_fragments,
            )
            for item in lookup(dpath="simulation/custom", within=config)
            for n_fragments in item["n_fragments"]
        ]

    random.seed(lookup(dpath="simulation/random", within=config)["seed"])
    if lookup(dpath="simulation/random/strata", within=config) is not None:
        retval += [
            collect(
                patterns,
                seq=generate_random_sequences(
                    seq_len=random.choice(range(10, 21)),
                    mutation_rate=lookup(dpath="simulation/random", within=config)[
                        "mutation_rate"
                    ],
                    modified_nucleosides=workflow.source_path(
                        "../resources/masses.tsv"
                    ),
                ),
                n_fragments=n_fragments,
            )
            for seq in lookup(dpath="simulation/random/strata", within=config)
            for n_fragments in seq["n_fragments"]
            for _ in range(seq["num_seq"])
        ]

    if len(retval) == 0:
        print("No simulation data given.")

    return retval


def collect_experiments(*patterns):
    retval = []

    if lookup(dpath="experiment", within=config) is not None:
        retval += [
            collect(
                patterns,
                seq=item["seq"],
                n_fragments=item["fragments"],
            )
            for item in lookup(dpath="experiment", within=config)
        ]

    if len(retval) == 0:
        print("No experimental data given.")

    return retval
