from .utils import Atom, Residue, ActiveSite
aa3 = "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split()


def compute_similarity(site_a, site_b):
    """
    Compute the similarity between two given ActiveSite instances.

    Input: two ActiveSite instances
    Output: the similarity between them (a floating point number)
    """
    if len(site_a) > len(site_a):
        site_a, site_b = site_b, site_a
    dists = range(len(site_a) + 1)

    for i in site_b.index:
        nDists = [i +1]
        for j in site_a.index:
            aa_a = site_a.columns[site_a.iloc[j] == 1][0]
            aa_b = site_b.columns[site_b.iloc[i] == 1][0]
            if aa_a == aa_b:
                nDists.append(dists[j])
            else:
                m = min((dists[j], dists[j+1], nDists[-1]))
                nDists.append(1 + m)
        dists = nDists
    similarity = dists[-1]
    return similarity


def calc_avg_site_length(sites):
    ss = []
    for site in sites:
        ss.append(len(site.residues))

    return [sum(ss) / len(sites), max(ss), min(ss)]


def generate_random_site(sites):
    lens = calc_avg_site_length(sites)
    num_res = np.random.randint(lens[2],lens[1])
    site = pd.DataFrame(0, index = range(num_res), columns = aa3)

    for pos in site.index:
        aa = np.random.randint(0,19)
        site.iloc[pos,aa] = 1

    return site


def compute_cluster_center_dumb(cluster_list, sites):
    total = pd.DataFrame(columns = aa3)
    for j in cluster_list:
        site = sites[sites==j].onehot
        for row in site.index:
            if len(total) <= row:
                total = total.append(site.iloc[row])
            else:
                total.iloc[row] += site.iloc[row]
    return total / len(cluster_list)

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
            row.append(compute_similarity(i.onehot,j.onehot))
        simMat.append(row)
    return pd.DataFrame(simMat)

def cluster_by_partitioning(active_sites):
    """
    Cluster a given set of ActiveSite instances using a partitioning method.

    Input: a list of ActiveSite instances
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
    K-means algorithm, referenced from pythonprogramming.net
    https://pythonprogramming.net/k-means-from-scratch-machine-learning-tutorial/
    """
    cls = k_mean()
    cls = cls.fit_clusters(active_sites)

    return cls.classifications.values()


def cluster_hierarchically(active_sites):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  #

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """

    # Fill in your code here!

    return []




class k_mean:
    """
    document k mean class
    """
    def __init__(self, k =3, min_diff = 1, max_iter = 10):


        self.k = k
        self.min_diff = min_diff
        self.max_iter = max_iter


    def fit_clusters(self, sites):
            self.centroids = {}

            for i in range(self.k):
                self.centroids[i] = generate_random_site(sites, aa3)
            print('created random centroids')

            for i in range(self.max_iter):
                print(i)
                self.classifications = {}

                for i in range(self.k):
                    self.classifications[i] = []

                for site in sites:
                    distances = [compute_similarity(site.onehot, self.centroids[c]) for c in self.centroids]
                    classification = distances.index(min(distances))
                    self.classifications[classification].append(site)

                prev_centroids = dict(self.centroids)

                for classification in self.classifications:
                    self.centroids[classification] = compute_cluster_center_dumb(self.classifications[classification], sites)


                optimized = True

                for c in self.centroids:
                    original_centroid = prev_centroids[c]
                    current_centroid = self.centroids[c]
                    if  compute_similarity(original_centroid, current_centroid) > self.min_diff:
                        print(compute_similarity(original_centroid, current_centroid))
                        optimized = False

                if optimized:
                    print('opt')
                    return
