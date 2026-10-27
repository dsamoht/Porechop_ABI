"""
Ab-Initio inference of ONT adapter main file

https://github.com/qbonenfant
https://github.com/bonsai-team

This script evolved from the wrapper used to call the approx_counter
algorithm, which end goal is to find adapter sequence k-mers.
We are able to infer adapter sequence from the reads instead of using
the adapter.py static database (unmaintained since 2018) of native Porechop.
In order to change default parameters for approx_counter, change the config
file porchop/abinitio.conf or pass a custom config file based on the same
syntax using -abc / --ab_initio_config

This works has been funded by the french National Ressearch Agency (ANR).

This file is part of Porechop_ABI. Porechop_ABI is free software:
you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either
version 3 of the License, or (at your option) any later version. Porechop_ABI
is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with
Porechop_ABI. If not, see <http://www.gnu.org/licenses/>.
Author: Quentin Bonenfant (quentin.bonenfant@gmail.com)
"""

import os
from multiprocessing import cpu_count
from .adapters import Adapter
from collections import defaultdict as dd
import networkx as nx
import sys
from .arg_parser import get_arguments
from .drop_cut import *
from .consensus import *
from .cpp_function_wrappers import approx_counter
from math import inf

##############################################################################
#                                 CONSTANTS                                  #
##############################################################################

METHODS = ["greedy", "heavy"]  # Reconstruction methods.
ENDS = ["start", "end"]        # Possible ends key.
LOW_FREQ_THRESHOLD = 10        # K-mer low frequency treshold.


##############################################################################
#                              UTILITY FUNCTIONS                             #
##############################################################################


def warning(msg):
    """
    Print a formatted message to stderr
    @param msg: the message to print as warning (str)
    """
    print(f"/!\\\t {msg}", file=sys.stderr)


def low_freq_warning(end, amount):
    warning(f"The most frequent kmer has been found in less" +
            f" than {LOW_FREQ_THRESHOLD}% of the reads {end}s" +
            f" after approximate count ({amount}%)")
    warning("It could mean this file is already trimmed or the sample do not" +
            " contains detectable adapters.")


def print_consensus_result(consensus_adapters, v, print_dest=sys.stderr):
    """Display the result from consensus runs. The format is different
    beacause methods are pooled in consensus runs, so a method oriented
    display no longer makes sens.
    @param an adapter dictionnary, using appropriate ends as keys
    @param v, the verbosity level
    @param print destination, can be custom, but stdout is default
    """

    if(v > 0):
        print("\nRebuild adapters:\n",
              file=print_dest)

    for end in ENDS:
        print(end.capitalize(), file=print_dest)
        # TODO: Change the sorting method...
        # Using adapter name is only asking for problems.
        sorted_consensus = sorted(consensus_adapters[end].keys(),
                                  key=lambda x: int(x.split("_")[1]))
        for name in sorted_consensus:
            adp = consensus_adapters[end][name]
            print(name, file=print_dest)
            print(adp, file=print_dest)


def print_result(adapters, v, print_dest=sys.stderr):
    """Display the result in the legacy format (method->ends)
    using the updated data structure (end -> methods)
    @param an adapter dictionnary, using appropriate methods as keys
    @param v, the verbosity level
    @param print destination, can be custom, but stdout is default
    """
    out = print_dest

    if(v > 0):
        print("\n\nINFERRED ADAPTERS:\n",
              file=print_dest)

    # Getting methods from the dict directly.
    methods = sorted(adapters[ENDS[0]].keys())

    for meth in methods:
        msg = meth
        srt = adapters["start"][meth]
        end = adapters["end"][meth]
        if(v > 0):
            meth = meth.capitalize()
            srt = "Start:\t" + srt
            end = "End:\t" + end
        print(msg, file=print_dest)
        print(srt, file=print_dest)
        print(end, file=print_dest)


def print_adapter_dict(adapters, print_dest=sys.stderr):
    """Display the content of adapter counting dictionnaries.
    using the updated dict format.
    @param an adapter dictionnary, using appropriate methods as keys
    @param print destination, can be custom, but stderr is default
    """

    methods = sorted(adapters[ENDS[0]].keys())

    out = print_dest
    for meth in methods:
        srt, end = "", ""
        msg = meth
        for k, v in sorted(adapters["start"][meth].items(),
                           key=lambda x: adapters["start"][meth][x],
                           reverse=True):
            srt += f"{k}: {v}\n"
        for k, v in sorted(adapters["end"][meth].items(),
                           key=lambda x: adapters["end"][meth][x],
                           reverse=True):
            end += f"{k}: {v}\n"

        meth += " assembly method"
        meth = meth.capitalize()
        srt = "Start:\n" + srt
        end = "End:\n" + end

        print(msg, file=print_dest)
        print(srt, file=print_dest)
        print(end, file=print_dest)


def parse_conf(conf_file):
    """
    Parse the configuration file passed to approxCounter.
    It allows Porechop_ABI to know which parameters are
    used during the approximate k-mer count.
    @param The path to the configuration file (str)
    @return A dictionnary containing the parameters.
    """

    # Setting some default values:
    params = {'k': 16,      # kmer size, 2<= k <= 32
              'sl': 100,     # sequence sampling size
              'sn': 40000,   # number of sequence sampled
              'limit': 500,  # number of kmers to keep.
              }
    try:
        assert(os.path.isfile(conf_file))
    except AssertionError:
        warning("Invalide configuration file.")
        warning(conf_file)
        sys.exit(1)
    else:
        with open(conf_file) as f:
            for line in f:
                line = line.rstrip("\n").replace(" ", "")
                data = line.split("=")
                if(line and line[0] != "#" and len(data) == 2):
                    p, v = data
                    # int values
                    if(v.isdigit()):
                        params[p] = int(v)
                    # float values
                    elif(v.replace(".","",1).isdigit()):
                        params[p] = float(v)
                    # string values
                    else:
                        params[p] = v
    return(params)


##############################################################################
#                           ADAPTER SEQUENCE ASSEMBLY                        #
##############################################################################

def greedy_assembly(g):
    """Greedy assembly method to compute the adapter using the graph as input.
    @param  g: A NetworkX DiGraph (graph)
    @return path: The longest extendable sequence from the most frequent kmer.
    """

    # Finding the best node in the graph.
    start = max(g.nodes, key=lambda x: g.nodes[x]["weight"])
    path = [start]

    right_node = start
    left_node = start
    while(left_node or right_node):

        # forward extension
        if(right_node):
            r_list = [el for el in g.successors(right_node) if el not in path]
            if(not r_list):
                right_node = None
            else:
                right_node = max(r_list, key=lambda x: g.nodes[x]["weight"])
                path.append(right_node)

        # reverse extension
        if(left_node):
            l_list = [el for el in g.predecessors(left_node) if el not in path]

            if(not l_list):
                left_node = None
            else:

                left_node = max(l_list, key=lambda x: g.nodes[x]["weight"])
                path = [left_node] + path
    return(path)


def dag_heaviest_path(G):
    """Returns the heaviest path in a DAG

    @param G: A NetworkX DiGraph (graph)
    @return path: the heaviest path in the graph (list)

    // Note: //
    This is a modified version of the dag_longest_path() from Networkx
    using node weight as distance instead of edge weight.
    _________________________________________________________________
    """
    dist = {}  # stores [node, distance] pair
    for node in nx.topological_sort(G):
        # pairs of dist,node for all incoming edges
        pairs = [(dist[v][0] + G.nodes[v]["weight"], v) for v in G.pred[node]]
        if pairs:
            dist[node] = max(pairs)
        else:
            dist[node] = (0, node)
    node, (length, _) = max(dist.items(), key=lambda x: x[1])
    path = []
    while length > 0:
        path.append(node)
        length, node = dist[node]
    return list(reversed(path))


def heavy_path(g):
    """ Searching the heaviest path between all source and target nodes.
        Even if longer path tend to be heavier, very heavy short path can
        also be selected.
        @param g: A NetworkX DiGraph (graph)
        @return hv_path : Heaviest path in the graph using node weight (list)
    """

    hv_path = dag_heaviest_path(g)

    # Annotating graph
    for n in hv_path:
        g.nodes[n]["path"] = (g.nodes[n]["path"] + ",heavy").lstrip(",")

    return(hv_path)


def greedy_path(g):
    """Greedy assembly of the adapter from the graph,
    starting from the most frequent
    @param  g: A NetworkX DiGraph (graph)
    @return path: The longest extendable sequence from the most frequent kmer.
    """

    gd_path = greedy_assembly(g)
    # annotating graph
    for n in gd_path:
        g.nodes[n]["path"] = (g.nodes[n]["path"] + ",greedy").lstrip(",")
    return(gd_path)


##############################################################################
#                              GRAPH FUNCTIONS                               #
##############################################################################

def haveOverlap(seq1, seq2):
    """Check if the sequence 1 is a k-1 prefix of sequence 2
    @param first sequence
    @param second sequence
    @return seq1 and seq2 merged in one sequence on the overlap.
    """
    minOverlap = min(len(seq1), len(seq2)) - 1
    if(seq1[-minOverlap:] == seq2[:minOverlap]):
        return(seq1 + seq2[minOverlap:])
    else:
        return("")


def get_weight(g, path):
    """Compute the weight of a path in an Nx graph
    @param g: The Networkx DiDraph (graph)
    @param path: The path (list)
    @return total: The total weight of this path (int)
    """
    total = 0
    for node in path:
        total += g.nodes[node]["weight"]
    return(total)


def concat_path(path):
    """Concat the kmers of a path into a single sequence
    @param The path as a k-mer list (list)
    @return the full sequence. (string)
    """
    return(path[0][:-1] + "".join(el[-1] for el in path))


def build_graph(count_file):
    """Build the adapter from the kmer count file.
       The way it is done is by building a directed weighted graph
       and searching for the heaviest path.
       I also added the greedy adapter output
       @param count_file : the path to the file to build the graph (string)
       @return g: A Networkx DiGraph containing all weighted k-mers (graph)
    """
    kmer_count = {}
    kmer_list = []

    g = nx.DiGraph()

    # Try to build the graph.
    try:
        with open(count_file, 'r') as f:
            for line in f:
                km, nb = line.rstrip("\n").split("\t")
                kmer_count[km] = int(nb)
                kmer_list.append(km)
                g.add_node(km, weight=int(nb))  # , path = "" )

    # If the file can not be opened.
    except FileNotFoundError:
        warning("Unable to open k-mer count file:")
        warning(count_file)
        warning("Either the file was moved, deleted, or filename is invalid.")

    # in case we try to convert to int something that is not a number
    except TypeError:
        warning("Something is wrong with the k-mer count file format:")
        warning(count_file)
        warning("Each line need to be in this format: <k-mer> \\t <count>")

    else:
        # searching overlaps
        for km in kmer_list:
            for km2 in kmer_list:
                # since i do all the possible combinations
                # There is no need to test both orientation
                direct = haveOverlap(km, km2)
                if(direct != ""):
                    g.add_edge(km, km2)

        # removing singletons
        g.remove_nodes_from([node for node in g.nodes if len(
            list(nx.all_neighbors(g, node))) == 0])

        # Returning only the biggest connected component
        g = g.subgraph(max(nx.weakly_connected_components(g),
                           key=lambda x: get_weight(g, x)))
    finally:
        return(g)

##############################################################################
#                              ADPATER BUILDING                              #
##############################################################################


def execApproxCounter(args, out_file, conf_file, mr, v, print_dest=sys.stderr):
    """ Execute the approximate k-mer counting (in-process, same output
    files as the approx_counter program).
    @param args: The arguments passed to Porechop (argparse arguments)
    @param out_file : The prefix for the k-mer count output file.
    @param mr : number of run to perform (for multi run). default is 1 (int)
    @param v: verbosity level. Default is 1. (int)
    @param print_dest : File to print to, std out by default (file stream)

    // Note: //
    This function returns nothing, it either generate a k-mer counting file
    or stops the whole program.
    ________________________________________________________________________
    """

    # TODO: Manage number of threads in a better way.
    nb_thread = min(args.threads, cpu_count())

    if(v > 0):
        print("Using config file:" + conf_file, file=print_dest)

    # Counting in-process (same output files as the approx_counter program)
    print_dest.flush()
    sys.stdout.flush()
    seed = getattr(args, "seed", None)
    if(seed is not None):
        # Each batch gets its own stream, so the "_sup" batch differs from the first one.
        seed += out_file.count("_sup")
    approx_counter(args.input, conf_file, out_file, mr, v, nb_thread, seed)


def generic_build(method_fct, name, g, adapters, which_end, v, w, print_dest,
                  drop_cut=True):
    """Generic function used to build adapter from a path in the graph g.
    This function need a path building function (which returns a list of
    consecutive overlapping k-mers). The resulting adapter is then
    added in the supplied adapters dictionnary for further processing.
    @param method_fct: Handle to the function used to build the path. (python function)
    @param name : The name of the method used (string)
    @param g :  A NetworkX DiGraph (graph)
    @param adapters: An adapter dictionnary to store results. (defaultdict)
    @param which_end: Name of the sequence end currently processed: start or end (string)
    @param v : verbosity level. Default is 1. (int)
    @param w : Window size for the drop cut mean smoothing
    @param print_dest : File to print to, std out by default (file stream)
    """
    if(v >= 1):
        print(f"\tBuilding {name} {which_end} adapter", file=print_dest)

    try:
        method_p = method_fct(g)
    except(nx.exception.NetworkXUnfeasible,
           nx.exception.NetworkXNotImplemented,
           nx.exception.HasACycle) as nxE:

        warning(f"Could not compute {which_end} adaper with {name} method")
        warning("The resulting graph probably contains a loop.")
        warning("For more details, read the error report below:")
        warning(nxE)
        warning("\n")
        adapters[which_end][name] = ""

    else:
        cut_method_p = []
        if(not drop_cut):
            cut_method_p = method_p
        elif(which_end == "start"):
            cut_method_p = start_cut(method_p, g, w)
        elif(which_end == "end"):
            cut_method_p = end_cut(method_p, g, w)
        # If the method path is not empty, concat it in a string
        if(cut_method_p):
            adapters[which_end][name] = concat_path(cut_method_p)
        else:
            warning(f"Something went wrong for {which_end} \
                  adaper with {name} method")
            warning("The returned graph path is empty")
            adapters[which_end][name] = ""


def make_adapter_objects(adapters, v, print_dest):
    adp = []
    # If we just need to print the adapter
    if(v >= 1):
        print("Building adapter object", file=print_dest)

    methods = sorted(adapters[ENDS[0]].keys())

    for meth in methods:
        adp.append(
            Adapter(f"abinitio_{meth}_adapter",
                    start_sequence=(f'abinitio_{meth}_Top',
                                    adapters['start'][meth]),
                    end_sequence=(f'abinitio_{meth}_Bottom',
                                  adapters['end'][meth])))

    if(v >= 1):
        print("The inference of adapters sequence is done.",
              file=print_dest)
    return(adp)


def build_adapter(args, out_filename, nb_reads, v, print_dest):
    """Building adapters from counts using different method:
        - Building a Debruijn graph and searching heaviest path
        - Greedy assembly based on kmer rank (most frequent first)
    Porechop can still properly detect adapters if there is a slight
    difference  since it uses a quite low identity threshold of 75%.
    @param A De Bruijn graph with kmer count as weight
    @param The argparse argument dictionnary from porechop.py
    @return A list of Adapter objects containing infered adapters.
    """

    # window size for the count smoothing at the drop cut stage.
    w = args.window_size

    # Keeping track of the best k-mer frequency
    best_freq = []

    # Placeholder adapter string dict
    adapters = dd(dict)  

    # Failsafe if the graphs for the adapter can not be built
    unable_to_build = False

    for which_end in ENDS:
        if(v >= 1):
            print("Assembling " + which_end + " adapters", file=print_dest)

        # Building graph
        g = build_graph(out_filename + "." + which_end)

        # Checking graph is properly build
        if(not nx.is_empty(g)):

            # Computing best k-mer frequency to decide
            # if a warning is needed or not
            max_count = max(g.nodes[n]["weight"] for n in g.nodes)
            best_freq.append(max_count)

            # preping for automatic anotation
            nx.set_node_attributes(g, "", "path")

            # Greedy Adapter
            generic_build(greedy_path, "greedy",
                          g,
                          adapters,
                          which_end,
                          v,
                          w,
                          print_dest,
                          not args.no_drop_cut)
            # Heavy adapter
            generic_build(heavy_path, "heavy",
                          g,
                          adapters,
                          which_end,
                          v,
                          w,
                          print_dest,
                          not args.no_drop_cut)

            # Exporting graph, if required
            if(args.export_graph is not None):
                if(v >= 1):
                    print("\tExporting assembly graph", file=print_dest)
                base = "_adapter_graph.graphml"
                if(".graphml" in args.export_graph):
                    base = "_" + os.path.basename(args.export_graph)
                path = os.path.join(os.path.dirname(
                    args.export_graph), which_end + base)

                nx.write_graphml(g, path)
        else:
            unable_to_build = True
            break

    if(unable_to_build):
        print("#################################", file=sys.stderr)
        print("ERROR  -  Unable to build graph  ", file=sys.stderr)
        print("Count file format may be invalid ", file=sys.stderr)
        print("#################################", file=sys.stderr)
        sys.exit(1)

    # Low frequency warnings
    lf_warning = [False, False]
    for i in range(len(ENDS)):
        current_freq = (best_freq[i] * 100.0) / nb_reads
        if(current_freq < LOW_FREQ_THRESHOLD):
            low_freq_warning(ENDS[i], current_freq)
            lf_warning[i] = True

    return (adapters, lf_warning)

##############################################################################
#                            Multi Run / Consensus                           #
##############################################################################


def insert_adapter_in_adp_count(adapter, adp_count):
    """Insert an adapter in an adapter collection dictionnary
    @param adapter A nested adapter dictionnary (start/end -> method used)
    @param adp_count A collection of previous adapter (including frequency)
    """
    for which_end in adapter.keys():
        methods = adapter[which_end].keys()
        for method in methods:    
            # adding ends dictionnary to adp_count
            if(which_end not in adp_count.keys()):
                adp_count[which_end] = {}
            # adding method sub_dictionnary to adp_count
            if(method not in adp_count[which_end].keys()):
                adp_count[which_end][method] = dd(int)

            # updating adapter count
            the_adapter = adapter[which_end][method]
            adp_count[which_end][method][the_adapter] += 1


def make_consensus_adapter_objects(adapters, v, print_dest):
    """ Build Adapter object with only start or end sequence.
    Pairing consensus sequences do not make sense since they
    come from different samples.
    @param adapters: The nested dictionnary with adapter sequences
    @param v: verbosity level (int)
    @param print_dest: the file to print the info (stream)
    """
    adp = []
    if(v >= 1):
        print("Building consensus adapter objects", file=print_dest)

    # start
    for consensus_name in adapters["start"].keys():
        adp.append(
            Adapter(f"{consensus_name}_adapter",
                    start_sequence=(f'{consensus_name}',
                                    adapters['start'][consensus_name])))
    # end
    for consensus_name in adapters["end"].keys():
        adp.append(
            Adapter(f"{consensus_name}_adapter",
                    end_sequence=(f'{consensus_name}',
                                  adapters['end'][consensus_name])))

    if(v >= 1):
        print("The inference of adapters sequence is done.",
              file=print_dest)
    return(adp)


def consensus_adapter(args, prefix, conf_file, nb_reads, v, print_dest):
    """ Try to find the fittest adapter for the dataset by running multiple
    adapter reconstruction and building a consensus.
    Default is to first try 10 runs, and keep the adapter if a perfectly stable
    consensus is found immediatly. Otherwise, 20 other runs are launched and
    the adapter sequence is build from a consensus of the 30 total runs.
    """

    # Final adapter object list
    adp = []
    # Adapters placeholders
    adapters = dd(dict)

    # Consensus adapter counter
    adp_count = {}

    # If we only need to print the inferred adapter
    just_print = args.guess_adapter_only

    # Defining output filename basename
    out_filename = prefix + "approx_kmer_count"

    #################################################################
    # CONSENSUS BUILDING
    #
    # Starting number of runs
    nb_run = args.number_of_run
    total_run = nb_run
    # While we have no consensus, build new runs.
    consensus_found = False
    # First batch need to have 100% consensus, else we rerun
    first_batch_done = False
    # Unless an error occured.
    build_error = False

    # Keeping track of the warnings
    triggered_warnings = {"LF": [0, 0], "CW": [0, 0]}

    # Changing verbosity level in multirun
    mr_v = v
    if(nb_run > 1 and v < 2):
        mr_v = 0
    while not consensus_found and not build_error:
        if(v > 0):
            if(not first_batch_done):
                print(f"Starting with a {nb_run} run batch.", file=print_dest)
            elif(nb_run > 0):
                print(f"Following with a {nb_run} run batch.", file=print_dest)

        # Launching and parsing results from approxCounter
        if(nb_run > 0):
            # Launching approx_counter in multi run mode
            execApproxCounter(args, out_filename, conf_file, nb_run, v, print_dest)

            # Parsing individual files
            for i in range(nb_run):
                print(f"Assembling run {i + 1}", file=print_dest)
                try:
                    current_adp, lf_warning = build_adapter(args,
                                                            out_filename + f"_{i}",
                                                            nb_reads,
                                                            mr_v,
                                                            print_dest)
                except FileNotFoundError:
                    # File not found...
                    warning("One of the k-mer count file was not found.")
                    warning("Something must have gone wrong with approxCounter.")
                    exit(1)

                # storing and counting adapters
                insert_adapter_in_adp_count(current_adp, adp_count)
                # If a low frequency warning is triggered, keep track of it
                triggered_warnings["LF"][0] += 1 if lf_warning[0] else 0
                triggered_warnings["LF"][1] += 1 if lf_warning[1] else 0

        # pooling methods for consensus mode.
        pooled_name = "consensus_pool"
        pooled = pool_methods(adp_count, pooled_name)
        # If this is the first run:
        if(not first_batch_done):
            first_batch_done = True
            consensus_found = True
            # Do we have a perfect consensus ?
            for end in ENDS:
                adp_list, counts = zip(*pooled[end][pooled_name].items())
                # if we have more than one adapter in the pooled dict, then no
                if(len(adp_list) != 1):
                    consensus_found = False
                else:
                    adapters[end][f"Adapter_1_{end}_(100%)"] = adp_list[0]
            if(consensus_found and v > 0):
                print("Perfect consensus found!\n"
                      f"Same sequence in {nb_run/nb_run} runs",
                      file=print_dest)

        # If we are at the second run, we need to find a consensus
        else:
            consensus_found = True
            # building the consensus dictionnary from the pooled adapter dict.
            adapters = build_consensus_adapter_dict(args,
                                                    pooled,
                                                    total_run * len(METHODS))
            # checking adapter list is empty
            if(len(adapters) == 0):
                warning("Adapter list is empty after first pass.")
                warning("Something very wrong happened and need to be fixed.")
                warning("If you were using the program 'as intended',\
                        please open an issue with details about this run\
                        on the Porechop_ABI github page.")
                exit(1)

        # Do we need to rerun ?
        if(not consensus_found):
            # If a 100% consensus is not found at first,
            # we increase the number of run and add those results
            # to current adapter collection (default: 20)
            nb_run = args.consensus_run
            total_run += nb_run
            out_filename += "_sup"
            warning("The adapters currently found are not all identical.")
            warning("A consensus will be made.")
            if(nb_run > 0):
                warning(f"Launching {nb_run} additional runs "
                        "to build a consensus.")

    # In case we want to export temp adapters used for the consenus
    if(args.export_consensus):
        try:
            out = open(args.export_consensus, "wt")
        except FileNotFoundError:
            warning("Could not export consensus file to")
            warning(args.export_consensus)
        else:
            print_adapter_dict(pooled, out)
            out.close()

    if(v > 0):
        print("Consensus step done", file=print_dest)

    #################################################################
    # RESULT EXPORT
    #
    # Are their any frequency warnings ?
    for i in range(len(ENDS)):
        nb_warn = triggered_warnings["LF"][i]
        if(nb_warn > 0):
            warning(f"Frequency warning triggered by {nb_warn}/{total_run} "
                    f"runs for {ENDS[i]} adapters")

    #########################################################
    # DEPRECATED: Old code just printed if -go was set.
    # Now printing adapters either way.
    #############
    # # If we just need to print the adapter
    # if(just_print):

    print_consensus_result(adapters, v, print_dest)

    # If we need to use them, we build porechop Adapter objects
    if(not just_print):
        adp = make_consensus_adapter_objects(adapters, v, print_dest)
    return adp


##############################################################################
#                                  MAIN                                      #
##############################################################################


def launch_ab_initio(args):
    """ Launch approx_counter to perform the approximate kmers count and try
    to rebuild the adapter using different methods.
        The count will be stored as a temporary file in a ./tmp folder.
        @param path to fasta file.
        @return the potential adapters as an Adapter object.
    """

    # getting args
    print_dest = args.print_dest
    v = args.verbosity

    # Preparing temporary file folder
    filename_pref = os.path.join(args.temp_dir, "temp_")
    tmpDir = os.path.dirname(filename_pref)

    # Creating tmp folder if it does not exist.
    if(not os.path.exists(tmpDir)):
        try:
            os.makedirs(tmpDir, exist_ok=True)
        except OSError as e:
            warning(e)
            warning("")
            warning("The tmp folder required for k-mer count files "
                    "could not be created.\n"
                    "Please check the error message above and provide "
                    "a writable path using the -tmp option.")
            exit(1)

    # Check if folder is writable:
    if(not os.access(tmpDir, os.W_OK)):
        warning("The tmp folder required for Porechop_ABI is not writable")
        warning("Please provide a writable path using the -tmp option.")
        exit(1)

    # temp adapter object placeholder
    adp = []

    # Using custom config file if specified
    conf_file = os.path.join(os.path.dirname(
        os.path.realpath(__file__)), "ab_initio.config")
    conf_file = args.ab_initio_config if args.ab_initio_config else conf_file

    # getting params sent to approxCounter
    params = parse_conf(conf_file)

    # Getting real number of read by comparing sample size and
    # number of reads in the fasta/q
    nb_reads = min(args.nb_reads, params["sn"])

    # Executing consensus adapter inference algorithm
    adp = consensus_adapter(args,
                            filename_pref,
                            conf_file,
                            nb_reads,
                            v,
                            print_dest)

    #########################################################
    # DEPRECATED: previously separated single /multirun mode.
    #############
    # else:
    #     # Or regular adapter inference algorithm
    #     out_filename = filename_pref + "approx_kmer_count"

    #     # Launching approx_counter in simple mode
    #     execApproxCounter(args, out_filename, conf_file, 1, v, print_dest)

    #     # Finally building  the adapter
    #     adapters, lf_warning = build_adapter(args, out_filename, nb_reads, v, print_dest)

    #     # If we need to just print the adapter, do it
    #     if(args.guess_adapter_only):
    #         print_result(adapters, v, sys.stdout)
    #     # else export adapter as Adapter object
    #     else:
    #         adp = make_adapter_objects(adapters, v, print_dest)
    #########################################################
    return(adp)


if __name__ == '__main__':
    args = get_arguments()
    adapt = launch_ab_initio(args)
    for ad in adapt:
        print(ad.__dict__)
