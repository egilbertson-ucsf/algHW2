from .utils import Atom, Residue, ActiveSite

def compute_similarity(site_a, site_b):
    """
    Compute the similarity between two given ActiveSite instances.

    Input: two ActiveSite instances
    Output: the similarity between them (a floating point number)
    """
    if len(site_a.residues) > len(site_b.residues):
        site_a, site_b = site_b, site_a
    dists = range(len(site_a.residues) + 1)

    for i, c in enumerate(site_b.residues):
        nDists = [i +1]
        for j, d in enumerate(site_a.residues):
            if c.type == d.type:
                nDists.append(dists[j])
            else:
                m = min((dists[j], dists[j+1], nDists[-1]))
                nDists.append(1 + m)
        dists = nDists
    similarity = dists[-1]
    return similarity

def calc_similarity_matrix(sites):
    """
    Calculate a complete matrix of similarities of all active sites to all others
        to be used to pull from in clustering so calculations only need to be done once
    Input: a list of ActiveSite instances
    Output: complete all by all matrix of levenstein distances between active sites
            formatted as a pandas DataFrame
            rows, columns = [0, 1, 2, ... n]

    """
    simMat = []
    for i in sites:
        row = []
        for j in sites:
            row.append(leven_dist(i,j))
        simMat.append(row)
    return pd.DataFrame(simMat)

def cluster_by_partitioning(active_sites):
    """
    Cluster a given set of ActiveSite instances using a partitioning method.

    Input: a list of ActiveSite instances
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
    """
    # Fill in your code here!

    return []


def cluster_hierarchically(active_sites):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  #

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """

    # Fill in your code here!

    return []
