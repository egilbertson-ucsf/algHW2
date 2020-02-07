from hw2skeleton import cluster as cl
from hw2skeleton import io
import sklearn.metrics as sk
import os
import pandas as pd
import numpy as np
import math
aa3 = "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split()
aa_df = pd.DataFrame(0, index=list(aa3), columns=['Count'])


def calc_avg_site_length(sites):
    '''
    calculate the average size of an active site
    for use in generating random sites
    '''
    ss = []
    for site in sites:
        ss.append(len(site.residues))

    return [sum(ss) / len(sites), max(ss), min(ss)]


def generate_random_site(sites):
    '''
    generate a random site by filling in a 1x20 vector repr of amino acids with counts
    '''
    lens = calc_avg_site_length(sites)
    num_res = np.random.randint(lens[2],lens[1])
    site = aa_df.copy()

    for pos in range(num_res):
        aa = np.random.randint(0,19)
        site.iloc[aa] += 1

    return site

def generate_k_random_centroids(k, sites):
    '''
    generate k random sites using above function
    '''
    centroids = {}
    for i in range(k):
        centroids[i] = generate_random_site(sites)
    return centroids

def assign_single_site_to_cluster(site, centroids):
    '''
    check which cluster centroid is closest to the given site and assign the
    site to that cluster
    '''
    loc = site.counts
    dists = {}
    for c in centroids.keys():
        dist = cl.compute_similarity(loc, centroids[c])
        dists[dist] = c
    closest = dists[min(dists.keys())]
    return closest

def assign_all_sites_to_cluster(sites, centroids, clusters):
    '''
    loop through all sites and assign them to the appropriate clusters
    '''
    for site in sites:
        close = assign_single_site_to_cluster(site, centroids)
        if close not in clusters:
            clusters[close] = [site]
        else:
            clusters[close].append(site)
    for cent in centroids:
        if cent not in clusters:
            clusters[cent] = []
    return clusters

def compute_cluster_center(cluster_list, sites_dict):
    '''
    compute the center of a cluster by taking the average of the vector representations
    of all sites in the cluster
    '''
    sites = aa_df.copy()
    for j in cluster_list:
        if isinstance(j, str):
            sites += sites_dict[j].counts
        else:
            sites += j.counts
    return sites / len(sites)

def get_new_centroids(clusters, sites_dict=None):
    '''
    use the compute_cluster_center function to get the new centroids after updating
    assignments
    '''
    centroids = {}
    for cluster in clusters.keys():
        centroids[cluster] = compute_cluster_center(clusters[cluster], sites_dict)
    return centroids

def check_change_in_centroids(old_centroids, new_centroids):
    ''' check how far the centroids have moved '''
    diff = 0
    for c in old_centroids.keys():
        diff += cl.compute_similarity(old_centroids[c], new_centroids[c])
    return diff

def one_full_k_means(sites, k):
    ''' using all above functions, one full iteration of k means'''
    centroids = generate_k_random_centroids(k, sites)
    clusters = {}
    clusters = assign_all_sites_to_cluster(sites, centroids, clusters)
    new_centroids = get_new_centroids(clusters)
    old_diff = check_change_in_centroids(centroids, new_centroids)
    new_diff = 0
    while old_diff - new_diff > 0.00001:
        old_diff = check_change_in_centroids(centroids, new_centroids)
        centroids = new_centroids.copy()
        clusters = {}
        clusters = assign_all_sites_to_cluster(sites, centroids, clusters)
        new_centroids = get_new_centroids(clusters)
        new_diff = check_change_in_centroids(centroids, new_centroids)
    return clusters, centroids

def compute_similarity_matrix(sites):
    ''' copy of computer similarity matrix from utils '''

    simMat = []
    names = []
    for i in range(len(sites)):
        names.append(sites[i].name)
        row = []
        for j in range(len(sites)):
            row.append(cl.compute_similarity(sites[i].counts,sites[j].counts))
        simMat.append(row)
    simMat = pd.DataFrame(simMat, columns = names, index = names)

    return simMat

def make_cluster_assign_df(clusters, simMat):
    ''' make a nice df repr of the cluster assignments'''
    assgn = pd.DataFrame(index = simMat.index, columns = ['Cluster Assignment'])
    for cluster in clusters.keys():
        for site in clusters[cluster]:
            assgn.loc[site.name] = cluster
    return assgn

def avg_sl(sites, k, simMat):
    ''' average silhouette_score for i random starts of k means for k clusters'''

    scores = []
    c_list = []
    for i in range(1):
        clusters, centroids = one_full_k_means(sites, k)
        assgn = make_cluster_assign_df(clusters, simMat)
        c_list.append(clusters)
        scores.append(sk.silhouette_score(simMat, assgn['Cluster Assignment'], metric='precomputed'))
    return scores, clusters



def k_means(sites=None):
    ''' run k means '''
    sites = io.read_active_sites('data')
    simMat = compute_similarity_matrix(sites)
    points = [[],[]]
    clusters = []
    for i in range(2,5):
        points[0].append(i)
        temp = avg_sl(sites, i , simMat)
        points[1].append(temp[0])
        clusters.append(temp[1])

    return clusters[points[1].index(max(points[1]))], max(points[1])
