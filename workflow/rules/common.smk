import re
import random

wildcard_constraints:
    n_fragments="[0-9]+",


_NUCLEOSIDE_RE = re.compile(r"\d*[ACGU]")


def get_seq_len(wildcards):
    return len(_NUCLEOSIDE_RE.findall(wildcards.seq))

def generate_random_sequences(seq_len,seed=0,mutation_rate=0):
    
    random.seed(seed)
    #Can add mutation rate here by defining the probabilities of the different bases (of course including the mutated bases)
    #Wait for the code for the modified bases!
    return ''.join([random.choice(["A", "U", "G", "C"]) for _ in range(seq_len)])


def collect_simulation_results(*patterns):
    return [
        collect(
            patterns,
            seq=item["seq"],
            n_fragments=n_fragments,
        )
        for item in lookup(dpath="simulation/custom", within=config)
        for n_fragments in item["n_fragments"]
    ] + [
    collect(
            patterns,
            seq=generate_random_sequences(
                seq_len = seq["seq_len"], 
                seed = lookup(dpath="simulation/random", within=config)["seed"], 
                mutation_rate = lookup(dpath="simulation/random", within=config)["mutation_rate"]
                ),
            n_fragments=n_fragments,
        )
        for seq in lookup(dpath="simulation/random/strata", within=config)
        for n_fragments in seq["n_fragments"]
        ]

