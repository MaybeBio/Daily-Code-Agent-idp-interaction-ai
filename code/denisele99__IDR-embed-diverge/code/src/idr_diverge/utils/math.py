import numpy as np


def transform_data(data_dict, type_):
    """
    Transforms a dictionary of lists based on the specified type.
    
    Supported types:
    - 'log': apply log10(x + 1) to each value
    - 'log_mean': apply log10(x + 1), then take the mean
    - 'geo_mean': apply log10(x + 1), then take 10 ** mean
    - 'mean': take the mean of each list
    """
    
    if type_ == 'log':
        return {k: np.log10([x + 1 for x in v]) for k, v in data_dict.items()}
    
    elif type_ == 'log_mean':
        return {k: np.mean(np.log10([x + 1 for x in v])) for k, v in data_dict.items()}
    
    elif type_ == 'geo_mean':
        return {k: 10 ** np.mean(np.log10([x + 1 for x in v])) for k, v in data_dict.items()}
    
    elif type_ == 'mean':
        return {k: np.mean(v) for k, v in data_dict.items()}
    
    elif type_ == 'median':
        return {k: np.median(v) for k, v in data_dict.items()}
    
    elif type_ == 'log_median':
        return {k: np.median(np.log10([x + 1 for x in v])) for k, v in data_dict.items()}
    
    elif type_ == 'geo_median':
        return {k: 10 ** np.median(np.log10([x + 1 for x in v])) for k, v in data_dict.items()}
    
    
    else:
        raise ValueError(f"Unsupported type: {type_}. Must be one of {'log','log_mean','geo_mean','mean', 'median','log_median','geo_median'}")
    

def permutation_differences(group_a, group_b, n_perm=10000):
    
    #group_a and group_b are lists of pairs/tuples
    
    # Combine data for permutation
    combined_data = np.concatenate((group_a, group_b))
    n_a = len(group_a)
    n_b = len(group_b)

    # 2. Calculate the observed test statistic (difference in means)
    observed_diff = np.mean(group_a) - np.mean(group_b)
    
    print(f"Observed difference in means: {observed_diff:.3f}")

    # 3. Perform permutations
    num_permutations = n_perm
    permutation_diffs = []

    for _ in range(num_permutations):
        # Shuffle the combined data
        np.random.shuffle(combined_data)
        
        # Split into new "group A" and "group B" based on original sizes
        perm_group_a = combined_data[:n_a]
        perm_group_b = combined_data[n_a:]    
        
        # Calculate the difference in means for this permutation
        perm_diff = np.mean(perm_group_a) - np.mean(perm_group_b)
        permutation_diffs.append(perm_diff)
    
    # 4. Calculate the p-value
    # Count how many permuted differences are as extreme as or more extreme than the observed difference
    # For a two-tailed test, we check both positive and negative extremes
    p_value = np.sum(np.abs(permutation_diffs) >= np.abs(observed_diff)) / num_permutations
    
    return permutation_diffs, p_value

def permutation_correlation(x:np.ndarray, y:np.ndarray, n_permutations: int = 10000):
    import numpy as np
    from scipy.stats import spearmanr

    # Calculate the observed Spearman correlation
    observed_corr, _ = spearmanr(x, y)
    
    count = 0

    # Perform permutation test by shuffling y
    for _ in range(n_permutations):
        y_permuted = np.random.permutation(y)
        perm_corr, _ = spearmanr(x, y_permuted)
        if abs(perm_corr) >= abs(observed_corr):
            count += 1

    # Calculate p-value
    p_value = count / n_permutations
    
    print(f"Observed Spearman correlation: {observed_corr:.4f}")
    print(f"Permutation test p-value: {p_value:.4f}")
    
    return observed_corr, p_value