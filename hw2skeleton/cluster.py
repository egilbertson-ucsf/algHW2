from .utils import Atom, Residue, ActiveSite
import numpy as np
import pandas as pd
aa3 = "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split()


def compute_similarity(site_a, site_b):
    """
    Compute the similarity between two given ActiveSite instances.

    Input: two ActiveSite instances
    Output: the similarity between them (a floating point number)
    """
    return np.linalg.norm(site_a - site_b)






def calc_similarity_matrix(sites):
    """
    Calculate a complete matrix of similarities of all active sites to all others
        to be used to pull from in clustering so calculations only need to be done once
    Input: a list of ActiveSite instances
    Output: complete all by all matrix of levenstein distances between active sites
            formatted as a pandas DataFrame
            rows, columns = [0, 1, 2, ... n]

    """

    return pd.DataFrame(simMat)

def cluster_by_partitioning(active_sites):
    """
    Cluster a given set of ActiveSite instances using a partitioning method.

    Input: a list of ActiveSite instances
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
    """


    return 


def cluster_hierarchically(active_sites):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  #

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """

    # Fill in your code here!

    return []



