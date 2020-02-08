from hw2skeleton import cluster as cl
from hw2skeleton import io
import os
import pandas as pd
import numpy as np
from hw2skeleton.utils import *
import sklearn.metrics as sk
aa3 = "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split()
aa_df = pd.DataFrame(0, index=list(aa3), columns=['Count'])


def compute_similarity_matrix(sites):
    '''
    Compute the similarity matrix of all sites by all sites using the similarity algorithm in the cluster.py script
    In: list of sites
    out: similarity matrix as pandas dataframe
    '''
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


def find_min(simMat):
    '''
    In: similarity matrix
    Out: a list of the names of the two clusters that are closest to one another
    simple loop through the similarity matrix
    '''
    min_value = simMat[simMat!=0].min()[0]
    min_pair = [simMat[simMat!=0].idxmin()[0], simMat.index[0]]
    for i in simMat:
        if simMat[simMat!=0][i].min() < min_value:
            min_pair = (simMat[simMat!=0][i].idxmin(), i)
            min_value = simMat[simMat!=0][i].min()
    return list(min_pair)


def rm_most_similar(df, pair):
    '''
    remove the rows and columns of the pair of most similar clusters from the similarity matrix

    '''
    return df.drop(pair, axis=0).drop(pair, axis =1)

def name_cluster(num):
    '''
    generate the string to name the nth clusters
    '''
    return str('c' + str(num))


def compute_cluster_center(cluster_list, sites_dict, aa_df):
    '''
    In list of sites in cluster, sites_dict to reference names to cluster objects, template DataFrame
    Out: vector repr of average of the cluster
    '''

    sites = aa_df
    for j in cluster_list:
        sites += sites_dict[j].counts
    return sites / len(sites)


def compute_new_cluster_sim(new_clust_avg, simMat_update, sites_dict, clusters):
    '''
    compute the similarity of the newly created cluster to the rest of the clusters
    '''
    newSim = []
    for site in simMat_update.columns:
        if site not in sites_dict:
            s = compute_cluster_center(clusters[site], sites_dict, aa_df)
        else:
            s = sites_dict[site].counts

        newSim.append(cl.compute_similarity(new_clust_avg, s))
    newSim.append(0.0)
    return newSim


def update_simMat(newSim, simMat_update, new_name):
    '''
    add the row/column of similarity values for the new cluster calculated above
    '''
    simMat_update[new_name] = None
    newRow = pd.DataFrame([newSim], columns = simMat_update.columns, index = [new_name])
    simMat_update = simMat_update.append(newRow)
    simMat_update[new_name] = newSim

    return simMat_update


def update_cluster_dict(new_name, min_pair, clusters, sites_dict):
    '''
    update the dictionary of cluster label to site lists
    '''
    pair0 = []
    pair1 = []
    if min_pair[0] in clusters:
        pair0 += unpack_cluster(clusters[min_pair[0]], sites_dict, clusters)
        del clusters[min_pair[0]]
    else:
        pair0.append(sites_dict[min_pair[0]].name)
    if min_pair[1] in clusters:
        pair1 += unpack_cluster(clusters[min_pair[1]], sites_dict, clusters)
        del clusters[min_pair[1]]
    else:
        pair1.append(sites_dict[min_pair[1]].name)
    clusters[new_name] = pair0 + pair1

    return clusters


def unpack_cluster(cluster_list, sites_dict, clusters):
    '''
    if a cluster name is not one of the base site names, unpack it to a list
    of the site names contained in the cluster
    '''
    out_list = []
    for key in cluster_list:
        if key not in sites_dict:
            for value in clusters[key]:
                out_list.append(value)
        else:
            out_list.append(key)
    return out_list

def run_everything(k, simMat_update, sites_dict):
    '''
    run all the above functions to get the cluster dictionary at a given level k
    '''
    c = 0
    clusters = {}
    while len(simMat_update) > k:
        new_name = name_cluster(c)
        c += 1
        min_pair=find_min(simMat_update)
        clust_sites = unpack_cluster(min_pair, sites_dict, clusters)
        simMat_update = rm_most_similar(simMat_update, min_pair)
        new_clust_avg = compute_cluster_center(clust_sites, sites_dict, aa_df)
        newSim = compute_new_cluster_sim(new_clust_avg, simMat_update, sites_dict, clusters)
        simMat_update = update_simMat(newSim, simMat_update, new_name)
        clusters = update_cluster_dict(new_name, min_pair, clusters, sites_dict)
    return clusters


def make_cluster_assign_df(clusters, simMat):
    '''
    helper function to make a nice representation for inputting into test_similarity
    '''
    assgn = pd.DataFrame(index = simMat.index, columns = ['Cluster Assignment'])
    for cluster in clusters.keys():
        for site in clusters[cluster]:
            assgn.loc[site] = cluster
    return assgn


def agglomerative(site = None):
    '''
    main function to run a complete agglomerative clustering job
    '''
    sites = io.read_active_sites('data')
    sites_dict = {}
    for site in sites:
        sites_dict[site.name] = site
    simMat = compute_similarity_matrix(sites)
    simMat_update = simMat.copy()
    c = run_everything(2, simMat_update, sites_dict)
    assgn = make_cluster_assign_df(c, simMat)
    return c, sk.silhouette_score(simMat, assgn['Cluster Assignment'], metric='precomputed')
