# Some utility classes to represent a PDB structure
import pandas as pd

## need this to create one hot encodings
aa3 = "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split()
def one_hot_encode(site):
    site.onehot = pd.DataFrame(0, index = range(len(site.residues)), columns=aa3)
    for i, aa in enumerate(site.residues):
        site.onehot.loc[i, aa.type] = 1
    return site.onehot

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
        self.onehot = None

    # Overload the __repr__ operator to make printing simpler.
    def __repr__(self):
        return self.name
