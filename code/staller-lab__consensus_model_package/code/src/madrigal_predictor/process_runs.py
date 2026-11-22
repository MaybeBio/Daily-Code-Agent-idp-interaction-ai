def get_runs(sequence, output, threshold):
    """Return contiguous runs of residues above and below threshold."""
    above_runs = []
    below_runs = []
    
    current_above = ''
    current_below = ''
    
    for aa, val in zip(sequence, output):
        if threshold is None or val >= float(threshold):
            current_above += aa
            if current_below:
                below_runs.append(current_below)
                current_below = ''
        else:
            current_below += aa
            if current_above:
                above_runs.append(current_above)
                current_above = ''
    
    # flush last run
    if current_above:
        above_runs.append(current_above)
    if current_below:
        below_runs.append(current_below)
    
    return above_runs, below_runs

def get_runs_coords(sequence, output, threshold):
    """Return contiguous runs of residues above and below threshold, and the coordinates of the transitions above and below the threshold"""

    above_runs = []
    below_runs = []
    starts_above = []
    starts_below = []

    current_above = ''
    current_below = ''

    for i in range(len(sequence)):
        if threshold is None or output[i] >= float(threshold):
            if not current_above:
                starts_above.append(i)
            current_above += sequence[i]
            if current_below:
                below_runs.append(current_below)
                current_below = ''

        else:
            if not current_below:
                starts_below.append(i)
            current_below += sequence[i]
            if current_above:
                above_runs.append(current_above)
                current_above = ''

    # flush last run
    if current_above:
        above_runs.append(current_above)
    if current_below:
        below_runs.append(current_below)


    return above_runs, below_runs, starts_above, starts_below
