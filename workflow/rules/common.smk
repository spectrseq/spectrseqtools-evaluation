import re


_NUCLEOSIDE_RE = re.compile(r"\d*[ACGU]")


def get_seq_len(wildcards):
    return len(_NUCLEOSIDE_RE.findall(wildcards.seq))


wildcard_constraints:
    n_fragments="[0-9]+",
