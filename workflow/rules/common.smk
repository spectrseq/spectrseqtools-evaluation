import re


wildcard_constraints:
    n_fragments="[0-9]+",


_NUCLEOSIDE_RE = re.compile(r"\d*[ACGU]")


def get_seq_len(wildcards):
    return len(_NUCLEOSIDE_RE.findall(wildcards.seq))


def collect_simulation_results(*patterns):
    return [
        collect(
            patterns,
            seq=item["seq"],
            n_fragments=n_fragments,
        )
        for item in lookup(dpath="simulation", within=config)
        for n_fragments in item["n_fragments"]
    ]
