'''
Get python modules
'''
from sys import argv
from math import sqrt, floor
from collections import Counter

'''
Get third-party modules
'''

'''
Get POSE modules
'''
from GetData import get_file

class PDB:

    """
    Reduce a pdb atom or hetatm line to a set of objects
    describing that atom
    """

    def __init__(self, Atom):
        
        self.AtomNumber = int(Atom[6:11]) 
        self.PDBAtom    = str(Atom[11:16]).strip()
        self.Element    = str(Atom[13:14])
        self.AminoAcid  = str(Atom[16:20]).strip()
        self.Chain      = str(Atom[20:22]).strip()
        self.Residue    = int(Atom[22:26])
        self.ChainRes   = str(Atom[20:22]).strip() + str(Atom[22:26]).strip()
        self.X          = float(Atom[30:38])
        self.Y          = float(Atom[38:46])
        self.Z          = float(Atom[46:54])
        self.XYZ        = [float(Atom[30:38]),float(Atom[38:46]),float(Atom[46:54])]

def center_of_geometry(Coordinates):
    '''
    Give me x, y, z coordinates, get back the cooresponding center of geometry
    '''
    
    return [sum(Coordinate)/len(Coordinate) for Coordinate in zip(*Coordinates)]

def distance(Coordinates):
    '''
    Euclidean distance between to x, y, z coordinates
    '''

    x,y,z = zip(*Coordinates)

    return sqrt((x[0] - x[1])**2 + (y[0] - y[1])**2 + (z[0] - z[1])**2)

def density(Distances, MinDist=0, MaxDist=10):
    '''
    For a given residues center of geometry, count the number of other residue
    center of geometries between a min and max disatance
    '''

    Counts = []
    
    for Distance in Distances:
       Counts.append(floor(Distance))
    
    Density = []
    for Distance, Count in sorted(Counter(Counts).items(), key=lambda Count: Count[0]):
        if Distance in range(MinDist, MaxDist):
            Density.append(Count)
                        
    return sum(Density)

def normalized_residue_burial(Arguments):
    '''
    This is just a very quick way to estimate which residues are buried and 
    and which are exposed. Basically get the atomic centers for 
    every residue (also could be a small molecule ligand) and count the 
    number of neighboring residues. Then normalize (i.e., the most buried
    residues will have a dad score of 1.0).
    '''

    #Get all PDB lines with atoms or heteroatoms
    Atoms = [Atom.strip() for Atom in file(get_file("Structure", Arguments.PathFile)) \
                 if Atom[:4] == "ATOM" or Atom[:6] == "HETATM"]
    #Make a list of residue (chain specific!) numbers and the corresponding aminoacids
    Residues, AminoAcids = zip(*sorted(set([(PDB(Atom).ChainRes, PDB(Atom).AminoAcid) for Atom in Atoms]), \
                                           key=lambda x: x[0]))

    #A dictionary of lists of lists, where every residues has all corresponding atomic coordinates
    Coordinates = dict([(Residue, []) for Residue in Residues])
    for Atom in Atoms:
        Coordinates[PDB(Atom).ChainRes].append(PDB(Atom).XYZ)

    #Center of geometry for every residue
    CenterOfGeometry = {}
    for Residue in Residues:
        CenterOfGeometry[Residue] = center_of_geometry(Coordinates[Residue])

    #For every residue, what's the distance to all other residues
    Distance = dict([(Residue, []) for Residue in Residues])
    for Residue1 in Residues:
        for Residue2 in Residues:
            Distance[Residue1].append(distance([CenterOfGeometry[Residue1], CenterOfGeometry[Residue2]]))

    #For every residue, how many residues are within some distance (here it's 0 to 10 angstroms)
    Density = dict([(Residue, []) for Residue in Residues])
    for Residue in Residues:
        Density[Residue] = density(Distance[Residue])

    AtomicDensityDict = {}
    Normalize = float(max(Density.values())) #Normalize by most densly packed residues
    for Residue in Residues:
        if Residue[0] == Arguments.Chain:
            AtomicDensityDict[int(Residue[1:])] = round(Density[Residue]/Normalize, 3)

    return AtomicDensityDict
        
            


                

