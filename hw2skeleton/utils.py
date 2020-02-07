# Some utility classes to represent a PDB structure
import pandas as pd

## need this to create one hot encodings
aa3 = "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split()
aa_df = pd.DataFrame(0, index=list(aa3), columns=['Count'])

def get_repr(site):
    df = aa_df.copy()
    for aa in site.residues:
        df.loc[aa.type] += 1
    return df


def make_repr_data_frame(sites, simMat):
    '''
    put the vector representations into a data frame for PCA
    '''
    df = pd.DataFrame(index = aa3, columns = simMat.columns )
    for site in sites:
        df[site.name] = site.counts
    return df

def make_cluster_assign_df(clusters, simMat):
    '''
    make a df of cluster assignments for downstream use
    '''
    assgn = pd.DataFrame(index = simMat.index, columns = ['Cluster Assignment'])
    for cluster in clusters.keys():
        for site in clusters[cluster]:
            assgn.loc[site.name] = cluster
    return assgn

def do_PCA(assgn, sites, simMat):
    '''
    Do PCA based on sklearn tutorial

    '''
    a = make_repr_data_frame(sites, simMat).T
    pca = PCA(n_components=2)
    p = pca.fit_transform(a)
    principalDf = pd.DataFrame(data = p
                 , columns = ['principal component 1', 'principal component 2'], index = simMat.index)
    finalDf = pd.concat([principalDf, assgn[['Cluster Assignment']]], axis = 1)

    return finalDF

def pca_plot(clusters, finalDf):
    '''
    plot pca based on sk learn tutorial
    '''
    fig = plt.figure(figsize = (8,8))
    ax = fig.add_subplot(1,1,1)
    ax.set_xlabel('Principal Component 1', fontsize = 15)
    ax.set_ylabel('Principal Component 2', fontsize = 15)
    ax.set_title('K-Means Clustering', fontsize = 20)
    targets = clusters.keys()
    colors = ['b', 'g']
    for target, color in zip(targets,colors):
        indicesToKeep = finalDf['Cluster Assignment'] == target
        ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1']
                   , finalDf.loc[indicesToKeep, 'principal component 2']
                   , c = color
                   , s = 50)
    ax.legend(targets)
    ax.grid()
    return

class Atom:
    """
    A simple class for an amino acid residue
    """

    def __init__(self, type):
        self.type = type
        self.coords = (0.0, 0.0, 0.0)

    # Overload the __repr__ operator to make printing simpler.
    def __repr__(self):
        return self.type

class Residue:
    """
    A simple class for an amino acid residue
    """

    def __init__(self, type, number):
        self.type = type
        self.number = number
        self.atoms = []

    # Overload the __repr__ operator to make printing simpler.
    def __repr__(self):
        return "{0} {1}".format(self.type, self.number)

class ActiveSite:
    """
    A simple class for an active site
    """

    def __init__(self, name):
        self.name = name
        self.residues = []
        self.counts = None

    # Overload the __repr__ operator to make printing simpler.
    def __repr__(self):
        return self.name
