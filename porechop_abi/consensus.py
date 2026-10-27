"""
Consensus building tools for Porechop_ABI

https://github.com/qbonenfant
https://github.com/bonsai-team

When using the multi run mode, several potential adapters are generated.
This file contains the required tools to build and manage the generated
sequences, clusterise them and output a consensus for each cluster of
compatible sequences.

This works has been funded by the french National Ressearch Agency (ANR).

This file is part of Porechop_ABI. Porechop_ABI is free software:
you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either
version 3 of the License, or (at your option) any later version. Porechop_ABI is
distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with
Porechop_ABI. If not, see <http://www.gnu.org/licenses/>.
Author: Quentin Bonenfant (quentin.bonenfant@gmail.com)
"""

import sys
import os
from collections import defaultdict as dd
import subprocess
from .cpp_function_wrappers import compatibility_matrix, msa_consensus


# Best consensus has to be above 30% or it will trigger a warning
BAD_CONSENSUS_THRESHOLD = 30


def run_subprocess(program, arguments):
    prog_path = os.path.join(os.path.dirname(
        os.path.realpath(__file__)), program)
    command_line = [prog_path] + arguments
    result = subprocess.run(command_line,
                           # capture_output=True,
                            encoding="utf-8",stderr=subprocess.PIPE,stdout=subprocess.PIPE)
    return(result.stdout)



def get_compatibles(mat):
    """ Find compatible subsets of sequences id from the similarity matrix
    @param mat: A similarity matrix (list)
    @return compat: A list of compatible sequences ids.
    """
    clustered = []
    compat = []
    nb_line = len(mat)
    i = 0
    while(len(clustered) < nb_line):
        if(i not in clustered):
            current = [i]
            clustered.append(i)
            for j in range(nb_line):
                if(i != j and j not in clustered):
                    if(abs(mat[i][j]) != 0):
                        # check if potental , new node is compatible.
                        if(est_clique(mat, current + [j])):                
                            clustered.append(j)
                            current.append(j)
            compat.append(current)
        i += 1
    return(compat)


def all_vs_all_matrix(sequences):
    """Build the distance matrix using the 'compatibility' algorithm."""
    return(compatibility_matrix(sequences))


def pool_methods(adapters, name):
    """
    Combine the adapters from different inferring methods into
    one global method.
    @param adapters An adapter dictionnary to concatenate.
    @param name The name of the new method
    """

    # New pooled adapter count (nested default dicts)
    pooled = dd(lambda: dd(lambda: dd(int)))
    for end in adapters.keys():
        for meth in adapters[end].keys():
            for adp, count in adapters[end][meth].items():
                pooled[end][name][adp] += count
    return(pooled)


def merge_adapter_set(adp1, adp2):
    """ Add the result of adp2 into adp1
    """
    for end in adp2.keys():
        for meth in adp2[end].keys():
            adp1[end][meth] = adp2[end][meth]


def seqan_msa_consensus(sequences, tmp_dir, end, adp_count):
    weighted = []
    methods = adp_count[end].keys()
    for seq in sequences:
        # get the number of time it should be repeated
        nb_time = 0
        for meth in methods:
            if(seq in adp_count[end][meth].keys()):
                nb_time += adp_count[end][meth][seq]
        weighted.append((seq, nb_time))
    # Same output as the msa_consensus program on the sequences repeated nb_time times
    return(msa_consensus(weighted))


def export_as_fasta(seq, filename):
    with open(filename, "w") as out:
        for i, s in enumerate(seq):
            out.write(f">Sequence_{i}\n{s}\n")


def print_mat(mat):
    for line in mat:
        print("".join([f"{el:>5}" for el in line]))


def est_clique(mat, node_set):
    missing = 0
    lns = len(node_set)
    origin = node_set[0]  # original node
    for i in range(lns):
        n1 = node_set[i]
        # # DEBUG
        # print(f"Testing node {n1}")
        for j in range(i+1, lns):
            n2 = node_set[j]
            # id there is a miss and n2 is not included in the main node
            if(mat[n1][n2] == 0 and mat[origin][n2] != 2):
                missing += 1
                # # DEBUG
                # print(f"Unlinked node: {n2}")
        if(missing >= 2):
            return(False)
    return(True)


def group_sequences(group, sequences):
    """Return the list of original sequences from a cluster list
    @param group: the list of sequences id that are in the same group (int)
    @param sequences: the list of all sequences for the current end. (list)
    @return The list of sequences corresponding to the group ids'. (list)
    """
    return(list(sequences[i] for i in group))


def group_weight(group, sequences, end, adp_count):
    """Compute the total weight of a sequence cluster
    using the total number of occurences in the original
    adapter count dictionnary.
    @param group: the list of sequences id that are in the same group (int)
    @param sequences: the list of all sequences for the current end. (list)
    @param end : The end we are working on (string)
    @param adp_count: A nested dictionnary containint adapter counts (dict)
    @return weight : The total weight of the cluster (int)
    """
    weight = 0
    for g in group:
        seq = sequences[end][g]
        for meth in adp_count[end].keys():
            if(seq in adp_count[end][meth].keys()):
                weight += adp_count[end][meth][seq]
    return(weight)


def adpcount_2_sequences(adp_count):
    """  Convert an adapter dictionnary to a sequences list,
    separated by ends. (methods are pooled)
    @param adp_count: A nested dictionnary containint adapter counts (dict)
    @return sequences: dictionnary of distincts sequences for both ends (dict)
    """
    pooled = pool_methods(adp_count, "pooled")
    sequences = dd(list)
    for end in pooled.keys():
        adp = pooled[end]["pooled"]
        # keeping the list sorted
        ordered_adp = sorted(list(adp.keys()),
                             key=lambda x: adp[x],
                             reverse=True)
        sequences[end] = ordered_adp
    return(sequences)


def sort_groups(groups, sequences, end, adp_count):
    """ Sort the groups by total weight, descending.
    @param group: the list of sequences id that are in the same group (int)
    @param sequences: the list of all sequences for the current end. (list)
    @param end : The end we are working on (string)
    @param adp_count: A nested dictionnary containint adapter counts (dict)
    @return sorted_g : The list of groups sorted by descending weithgt (list)
    """
    sorted_g = sorted(groups,
                      key=lambda x: group_weight(x, sequences, end, adp_count),
                      reverse=True)
    return(sorted_g)


def build_consensus_adapter_dict(args, adp_count, total_seq):
    """ Main function of this script. It clusterise and compute consensus
    for all sequences on both ends.
    @param adp_count: A nested dictionnary containing adapter counts (dict)
    @return consensus_count: a nested dictionnary with consensus for both ends
    """
    consensus_count = dd(lambda: dd(lambda: dd(int)))
    sequences = adpcount_2_sequences(adp_count)
    for end in adp_count.keys():
        # Alignement matrix and compatibility clustering.
        mat = all_vs_all_matrix(sequences[end])
        compatibles = get_compatibles(mat)

        # Sorting to get only the best consensus if needed.
        sorted_compat = sort_groups(compatibles, sequences, end, adp_count)

        # temp_dir to store fasta file
        tmp_dir = args.temp_dir

        # If the user only want the x best consensus
        box = args.best_of_x
        # Otherwise, get all consensus.
        nb_consensus = box if box != 0 else len(sorted_compat)
        # Computing best consensus.
        for i, group in enumerate(sorted_compat[:nb_consensus]):
            # fetching sequences and counts for this group
            group_seq = group_sequences(group, sequences[end])
            weight = group_weight(group, sequences, end, adp_count)
            consensus = seqan_msa_consensus(group_seq, tmp_dir, end, adp_count)
            # storing consensus in similar structure as adp_count,
            # for back-compatibility.
            frequency = round((100*weight) / total_seq, 1)

            # checking if best consensus is above 30% frequency
            if(i == 0 and frequency < BAD_CONSENSUS_THRESHOLD):
                print(f"/!\\\tThe best {end} consensus use less than 30%",
                      "of the inferred adapters set.",
                      file=sys.stderr)
                print("/!\\\tThe input file either do not contains adapters",
                      file=sys.stderr)
                print("/!\\\tor has a highly variable region (ex: barcode).",
                      file=sys.stderr)
                print("/!\\\tYou may want to re-run several time in -go mode",
                      file=sys.stderr)
                print("/!\\\tand input manually curated adapters with -cap.",
                      file=sys.stderr)

            # Only export consensus with high enough frequency (default is all)
            if(frequency >= args.all_above_x):
                # Relative frequency stored directly into the name
                consensus_count[end][f"Consensus_{i+1}_{end}_({frequency}%)"] = consensus
    return(consensus_count)
