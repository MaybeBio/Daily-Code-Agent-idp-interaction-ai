import numpy as np
from itertools import combinations
from pygosemsim import similarity

from idr_diverge.go_enrichment.go_enrichment_ import load_go_sim_graph


def jaccard_similarity(a,b): #jaccard similarity
    '''Jaccard similarity, where the input is a list of go terms
        Args:
            a(list): list of GO IDs
            b(list): list of GO IDs'''

    if len(a)==0 or len(b)==0:
        return np.nan

    intersection = len(set(a).intersection(b)) #overlapping non-zero positions
    union = len(set(a).union(b))
    #union = (len(a)+len(b)) - intersection #unique non-zero positions in A and B
 
    return intersection / union

def fold_enrichment(observed_ratio, background_n = 19459, NN=50):
    observed_ratio = tuple(observed_ratio.split('/'))
    n_term_observed,n_term_background = observed_ratio

    observed = int(n_term_observed)/NN
    expected = int(n_term_background)/background_n
    
    return round(observed/expected,1)


def wang_similarity(term1, term2,G_graph):
    return similarity.wang(G_graph, term1, term2)

G = load_go_sim_graph()

def go_bma_similarity(go_terms1, go_terms2):
    
    if (len(go_terms1)<1) or (len(go_terms2)<1):
        return np.nan
     
    # Calculate similarity using BMA method
    similarity_scores = []
    for term1 in go_terms1:
        try:
            best_match_score = max(wang_similarity(term1, term2, G) for term2 in go_terms2)
            similarity_scores.append(best_match_score)
        except:
            continue

    for term2 in go_terms2:
        try:
            best_match_score = max(wang_similarity(term2, term1, G) for term1 in go_terms1)
            similarity_scores.append(best_match_score)
        except:
            continue
    
    bma_similarity = sum(similarity_scores) / (len(go_terms1) + len(go_terms2))
    
    return bma_similarity

def go_bma_mean_similarity(go_terms_list):
    '''
    BMA mean similarity across list of go term lists (all pairs vs all pairs)
    
    :param go_terms_list: list of lists of go terms
    :return: mean BMA similarity
    '''
    #go_terms_list = list(set([term for sublist in go_terms_list for term in sublist]))
    if len(go_terms_list) <2:
        return np.nan
    else:
        sims = []
        for terms1,terms2 in combinations(go_terms_list, 2):
            #print(terms1, terms2)
            if len(terms1)>0 and len(terms2)>0:
                sims.append(go_bma_similarity(terms1, terms2))
        #print(sims)
        mean_similarity = np.mean(sims)
        return mean_similarity


def go_jaccard_similarity(go_terms1, go_terms2):
    #Similarity of two sets of GO terms
     # Convert lists to sets
    set1 = set(go_terms1)
    set2 = set(go_terms2)
    
    # Calculate the intersection and union of the sets
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    # Calculate Jaccard similarity
    if len(go_terms1) == 0 or len(go_terms2) == 0:
        return 0
    return intersection / union

def go_jaccard_divergence(go_terms1, go_terms2):
    return go_jaccard_similarity(go_terms1, go_terms2) - 1