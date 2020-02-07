from .utils import Atom, Residue, ActiveSite
import numpy as np
import pandas as pd
from sklearn.metrics import jaccard_similarity_score
from .k_means import *
from .agglomerative import *
aa3 = "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split()


def compute_similarity(site_a, site_b):
    """
    Compute the similarity between two given ActiveSite instances.

    This will calculate the Euclidean distance between the numeric vector
    representation of two active sites.
    The numeric vector representation is based on counts of each
    amino acid present at a given active site.


    Input: two ActiveSite instances
    Output: the similarity between them (a floating point number)
    """
    return np.linalg.norm(site_a - site_b)


def cluster_by_partitioning(active_sites):
    """
    Cluster a given set of ActiveSite instances using a partitioning method.
    Will call to k means clustering implementation that is housed in another script
    Input: a list of ActiveSite instances
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
    """
    cls, sc = k_means(active_sites)

    return cls


def cluster_hierarchically(active_sites):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  #

    Calls to agglomerative clustering algorithm housed in another script
    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """


    cls, sc = agglomerative(active_sites)

    return cls


def jacaard(clusters_a, clusters_b):
        """
        Jacaard score from sklearn
        Input: two lists of cluster assignments
        Output: score
                (each clustering is a list of lists of Sequence objects)
        """
    return jaccard_similarity_score(clusters_a, clusters_b)
