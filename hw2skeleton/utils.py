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
