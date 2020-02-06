from hw2skeleton import cluster as cl
from hw2skeleton import io
import sklearn.metrics as sk
import os
import pandas as pd
import numpy as np
aa3 = "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split()
aa_df = pd.DataFrame(0, index=list(aa3), columns=['Count'])


def calc_avg_site_length(sites):
    ss = []
    for site in sites:
        ss.append(len(site.residues))

    return [sum(ss) / len(sites), max(ss), min(ss)]


def generate_random_site(sites):
    lens = calc_avg_site_length(sites)
    num_res = np.random.randint(lens[2],lens[1])
    site = aa_df.copy()

    for pos in range(num_res):
        aa = np.random.randint(0,19)
        site.iloc[aa] += 1

    return site

def generate_k_random_centroids(k, sites):
    centroids = {}
    for i in range(k):
        centroids[i] = generate_random_site(sites)
    return centroids

def assign_single_site_to_cluster(site, centroids):
    loc = site.counts
    dists = {}
    for c in centroids.keys():
        dist = cl.compute_similarity(loc, centroids[c])
        dists[dist] = c
    closest = dists[min(dists.keys())]
    return closest

def assign_all_sites_to_cluster(sites, centroids, clusters):
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
    sites = aa_df.copy()
    for j in cluster_list:
        if isinstance(j, str):
            sites += sites_dict[j].counts
        else:
            sites += j.counts
    return sites / len(sites)

def get_new_centroids(clusters, sites_dict=None):
    centroids = {}
    for cluster in clusters.keys():
        centroids[cluster] = compute_cluster_center(clusters[cluster], sites_dict)
    return centroids

def check_change_in_centroids(old_centroids, new_centroids):
    diff = 0
    for c in old_centroids.keys():
        diff += cl.compute_similarity(old_centroids[c], new_centroids[c])
    return diff

def one_full_k_means(sites, k):
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
    assgn = pd.DataFrame(index = simMat.index, columns = ['Cluster Assignment'])
    for cluster in clusters.keys():
        for site in clusters[cluster]:
            assgn.loc[site.name] = cluster
    return assgn

def avg_sl(sites, k, simMat):
    scores = []
    for i in range(10):
        clusters, centroids = one_full_k_means(sites, k)
        assgn = make_cluster_assign_df(clusters, simMat)
        scores.append(sk.silhouette_score(simMat, assgn['Cluster Assignment'], metric='precomputed'))
    return min(scores)
