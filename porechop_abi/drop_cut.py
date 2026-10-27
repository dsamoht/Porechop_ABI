from statistics import median, mean

CUT_RATIO = 0.075  # How much random variation is allowed for drop cut.

def start_cut(start_adp_data, g, w=7):
    """ Try to find the best place to cut the raw adapter path
     to get an appropriate adapter.
    This function cut the end of forward adapters.
    @param The path, as a kmer list
    @param The De Bruijn graph.
    @param (opt) w, window size for sliding median smoothing (default: 7)
    @return The adjusted adapter path
    """

    # Data for Start Adapters
    start_dist = [g.nodes[el]["weight"] for el in start_adp_data]

    # computing start epsilon based on counts
    start_epsilon = CUT_RATIO * max(start_dist)

    # Sliding median method
    start_sl_md = sliding_median(start_dist, w)
    # # Sliding mean method
    # start_sl_md = sliding_mean(start_dist, w)

    # Simple derivative
    start_sl_md = simple_deriv(start_sl_md)

    # Finding cutting point
    # start median zone method
    smz = find_median_zone(start_sl_md, start_epsilon)

    # Finding the position to cut
    s_cut = len(start_adp_data) - 1
    if(smz):
        s_cut = min(max(smz), s_cut)

    # Returning the cut path
    return(start_adp_data[0:s_cut])


def end_cut(end_adp_data, g, w=7):
    """ Try to find the best place to cut the raw adapter path
     to get an appropriate adapter.
    This function cut the end of end adapters.
    @param The path, as a kmer list
    @param The De Bruijn graph.
    @param (opt) w, window size for sliding median smoothing (default: 7)
    @return The adjusted adapter path
    """

    # Fetching weight distribution
    end_dist = [g.nodes[el]["weight"] for el in end_adp_data]

    # Computing end epsilon based on counts
    end_epsilon = CUT_RATIO * max(end_dist)

    # Sliding median method
    end_sl_md = sliding_median(end_dist, w)
    # # Sliding mean method
    # end_sl_md = sliding_mean(end_dist, w)

    # Simple derivative
    end_sl_md = simple_deriv(list(reversed(end_sl_md)))

    # Finding cutting point:
    # end median zone method
    emz = find_median_zone(end_sl_md, end_epsilon)

    # re-reversing the end array
    end_sl_md.reverse()
    emz = [len(end_sl_md) - (s + 1) for s in emz]

    # Finding the position to cut
    e_cut = 0
    if(emz):
        e_cut = max(0, min(emz) + 1)

    # Returning the cut path
    return(end_adp_data[e_cut:])


def sliding_median(counts, w=7):
    """ Compute the median of a of each count value using
    a sliding window approach.
    @param counts a list of integer
    @param w the window size
    @return the median values, in a list.
    """
    lc = len(counts)
    offset = w // 2  # offset to work on "middle" position
    results = []
    for i in range(lc):
        start = max(0, i - offset)
        end = min(i + offset, lc)
        # getting on middle position
        med_rng = counts[start: end]
        # An empty window (window size 1) means no smoothing.
        results.append(median(med_rng) if med_rng else counts[i])
    return(results)


def sliding_mean(counts, w=7):
    """ Compute the mean of a of each count value using
    a sliding window approach.
    @param counts a list of integer
    @param w the window size
    @return the mean values, in a list.
    """
    lc = len(counts)
    offset = w // 2  # offset to work on "middle" position
    results = []
    for i in range(lc):
        start = max(0, i - offset)
        end = min(i + offset, lc)
        # getting on middle position
        med_rng = counts[start: end]
        results.append(mean(med_rng))
    return(results)


def simple_deriv(counts):
    """ Compute the difference between consecutives values in a list
    @param counts a list of integer
    @return The difference between a value and it's neighbourg.
    """
    lc = len(counts)
    results = [0]
    for i in range(1, lc, 1):
        val = counts[i] - counts[i - 1]
        results.append(val)
    return(results)


def find_median_zone(counts, epsilon):
    """ Find the position in the counts list where the
    counts raise above a threshold.
    @param counts the count list
    @param epsilon margin of acceptable error
    """
    c_median = median(counts)
    s_rng = 0  # start of low zone
    is_zone = False
    s_zone = []
    for i in range(len(counts) - 1, -1, -1):
        c = counts[i]
        if(c < c_median - epsilon):
            is_zone = True
            s_rng = i
        else:
            if(is_zone):
                is_zone = False
                s_zone.append(s_rng)
    return(s_zone)