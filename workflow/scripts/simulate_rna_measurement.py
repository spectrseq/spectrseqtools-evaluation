import sys

sys.stderr = open(snakemake.log[0], "w")

import polars as pl
import random, re
from scipy.stats import norm

# regex for separating given sequence into nucleosides
nucleoside_re = re.compile(r"\d*[ACGU]")

# input
true_sequence = nucleoside_re.findall(snakemake.wildcards.seq)
n_fragments = int(snakemake.wildcards.n_fragments)
nucleoside_masses = pl.read_csv(snakemake.input[0], separator="\t")
max_n_parts = int(snakemake.config["parameters"][0]["max_n_parts"])
rel_error_rate = snakemake.config["parameters"][0]["rel_error_rate"]
backbone_masses = snakemake.config["backbone_weights"][0]

# helper functions
#def random_fragment():
#    l = random.randint(0, len(true_sequence) - 1)
#    r = random.randint(l + 1, len(true_sequence))
#    return l, r

# Divide the sequence into a maximum given number of "max_n_parts". 
# At high energies the sequence can be divided into multiple parts! Sample uniformly how many parts it gets broken into.
def random_fragment(num_parts = max_n_parts):
    if num_parts == 1: 
        #If it breaks into one part, the entire sequence remains intact.
        breakagepoints = [int(0)]
    elif num_parts <= len(true_sequence): 
        #If the number of parts is less than or equal to the sequence length, then the sequence gets broken into the number of parts.
        breakagepoints = sorted(random.sample(range(1, len(true_sequence)), num_parts - 1))
    elif num_parts > len(true_sequence): 
        #The number of parts cannot be greater than the sequence length.
        raise ValueError("The number of parts cannot be greater than the sequence length!")
    return breakagepoints
# TODO: Implement that in some cases there is no base pair generated, but only the backbone with sugar etc?

#Select how many parts the sequence gets broken into. #It can happen that the entire sequence remains intact.
breakagepoints = [random_fragment(num_parts = random.randint(1,max_n_parts)) for _ in range(n_fragments)]

# Generate "left" and "right" neuclotides based on the breakagepoints, 
# which are the exact indices of the nucleotides in the generated fragments:
# NOTE: The change here from before: now left and right are the exact indices of the python list of nucleotides,
# with both left and right included in the fragment. If left = right, then its a single nucleotide.
def generate_left_right(breakagepoints):
    left,right = [],[]
    for exp in breakagepoints:
        left.append(0)
        for part in exp:
            if part != 0:
                right.append(part-1)
                left.append(part)
        right.append(len(true_sequence)-1)
    return left,right

l_breakage,r_breakage = generate_left_right(breakagepoints) 

# sample random fragments from true sequence
fragments = pl.from_records(
    list(zip(l_breakage,r_breakage)), 
    schema=["left", "right"],
    orient="row",
)

# TODO: Implement breaking at other points
#Middle -> -H + PO2
#5'end  -> correct
#3'end  -> +PO2
backbone_middle_addition = backbone_masses["phosphorus"] + 2*backbone_masses["oxygen"] - backbone_masses["hydrogen"]
backbone_start_addition  = 0.
backbone_end_addition    = backbone_masses["phosphorus"] + 2*backbone_masses["oxygen"]

# calculate true masses based on known nucleoside masses
true_fragment_masses = [
    sum(
        nucleoside_masses.filter(pl.col("nucleoside") == b)
        .select(pl.col("monoisotopic_mass"))
        .item()
        for b in true_sequence[fragment["left"] : fragment["right"] + 1]
    )
    for fragment in fragments.iter_rows(named=True)
]

#Copying the fragment masses to a new list to add the backbone masses to the fragments.
true_fragment_masses_with_backbone = [elem for elem in true_fragment_masses]

# Add backbone masses to the fragments, based on the position of the fragment in the sequence!
for iter,fragment in enumerate(fragments.iter_rows(named=True)):
    
    if fragment["left"] == 0:
        true_fragment_masses_with_backbone[iter] += backbone_start_addition
    elif fragment["left"] != len(true_sequence)-1:
        true_fragment_masses_with_backbone[iter] += backbone_middle_addition
    
    if fragment["right"] == len(true_sequence)-1:
        true_fragment_masses_with_backbone[iter] += backbone_end_addition
    elif fragment["right"] != 0:
        true_fragment_masses_with_backbone[iter] += backbone_middle_addition
    
    true_fragment_masses_with_backbone[iter] += (
        backbone_middle_addition*max(fragment["right"]-fragment["left"]+1-2,0))
    
    
# simulate observed masses with noise -- we use the relative error rate to simulate the noise!
fragment_noise = norm.rvs(scale=rel_error_rate, size=len(true_fragment_masses_with_backbone))
observed_fragment_masses = [
    max(mass + mass*noise, 0.0)
    for noise, mass in zip(fragment_noise, true_fragment_masses_with_backbone)
]

# compile final dataframe
fragments = fragments.with_columns(
    (pl.col("left") == 0).alias("is_start (5')"),
    (pl.col("right") == (len(true_sequence))-1).alias("is_end (3')"),
    (pl.col("right") == pl.col("left")).alias("single_nucleotide"),
    pl.Series(true_fragment_masses).alias("true_nucleoside_mass"),
    pl.Series(true_fragment_masses_with_backbone).alias("true_mass_with_backbone"),
    pl.Series(observed_fragment_masses).alias("observed_mass"),
)

# write simulation output
fragments.write_csv(snakemake.output[0], separator="\t")
