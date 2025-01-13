import re, random

wildcard_constraints:
    n_fragments="[0-9]+",

_NUCLEOSIDE_RE = re.compile(r"\d*[ACGU]")

def get_seq_len(wildcards):
    return len(_NUCLEOSIDE_RE.findall(wildcards.seq))

def generate_random_sequences(seq_len,seed=0,mutation_rate=0,modified_nucleosides=None):
    
    random.seed(seed)

    if modified_nucleosides is not None:
        
        with open(modified_nucleosides, 'r') as file:
            lines = file.readlines()[1:]
        modified_nucleoside_names = [line.split("\t")[0] for line in lines if line.split("\t")[0] not in ["U", "A", "G", "C"]]
        
        # The weights defined using the mutation rate  by defining the probabilities of the different bases (of course including the normal bases)
        weights = [(1.-mutation_rate/len(modified_nucleoside_names))/4]*4 + [mutation_rate/len(modified_nucleoside_names)]*len(modified_nucleoside_names)
        
        return ''.join(random.choices(["A", "U", "G", "C"] + modified_nucleoside_names, weights = weights, k = seq_len))
        
    else:
        return ''.join([random.choice(["A", "U", "G", "C"]) for _ in range(seq_len)])

def collect_simulation_results(*patterns):

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

    if lookup(dpath="simulation/random/strata", within=config) is not None:
        retval += [ 
            collect(
            patterns,
            seq=generate_random_sequences(
                seq_len = seq["seq_len"], 
                seed = lookup(dpath="simulation/random", within=config)["seed"], 
                mutation_rate = lookup(dpath="simulation/random", within=config)["mutation_rate"],
                #modified_nucleosides = workflow.source_path("../resources/masses_purines_modification.tsv")
                modified_nucleosides = workflow.source_path("../resources/masses_all.tsv")
                ),
            n_fragments=n_fragments,
        )
        for seq in lookup(dpath="simulation/random/strata", within=config)
        for n_fragments in seq["n_fragments"]
        ]
    
    if len(retval) > 0:
        return retval
    else:
        raise IndexError("Define at least one custom sequences in the config file or define the seq_len for random sequence generation")

