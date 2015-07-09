'''
Get python modules
'''
from os import environ

'''
Get third-party modules
'''

'''
Get POSE modules
'''

def get_file(DataType, PathFile):
    '''
    For a particular data type (e.g., mutation or sequence) find the 
    file containing said data using moca paths file
    '''

    Path = [Line.split()[2] for Line in open(PathFile).readlines() \
                if Line.strip() and Line.split()[0] == DataType][0]
    if Path[:7] == "${HOME}":
        Path = Path.replace("${HOME}", environ["HOME"])
    
    return Path

def get_sequences(Arguments):
    '''
    Takes a file of sequences in fasta format and makes a dictionary with gene keys and sequence values. 
    The sequence values are in the required "reduced fasta" format. This entails removing all gaps from 
    the reference sequence and removing the corresponding positions from the remaining sequences. The 
    reference gene is the gene under study (i.e., the one you are scoring mutations in). This gene/sequence 
    MUST be the first in your fasta file!!!
    '''

    Genes = [Gene.strip() for Gene in file(get_file("Fasta", Arguments.PathFile)) if Gene[0] == ">"]
    for Gene in Genes:
        vars()[Gene] = []

    for Line in file(get_file("Fasta", Arguments.PathFile)):
        if Line[0] == ">":
            Gene = Line.strip()
        else:
            vars()[Gene].extend(Line.strip())

    Alignment = []
    for Gene in Genes:
        Alignment.append(vars()[Gene])

    ReducedAlignment = []
    for Position in zip(*Alignment):
        if Position[0] != "-":
            ReducedAlignment.append(Position)

    ReducedAlignment = zip(*ReducedAlignment)

    Sequences = {}
    for Gene in Genes:
        Sequences[Gene] = "".join(ReducedAlignment[Genes.index(Gene)])
       
    return Sequences

def get_reference_gene(Arguments):
    '''
    The reference gene is the gene under study (i.e., the one you are scoring mutations
    in). This gene/sequence MUST be the first in your fasta file!!!
    '''

    return open(get_file("Fasta", Arguments.PathFile)).readlines()[0].strip()

def get_identities(ReferenceGene, Sequences, Arguments):
    '''
    For a given reference gene, this function calculates the sequence identity of all 
    "Sequences" relative to the reference. Returns a dictionary with gene keys and fractional
    identities for values
    '''
    
    Identities = {}
    if Arguments.IdentitySegment:
        UpStream = int(Arguments.IdentitySegment[0])
        DownStream = int(Arguments.IdentitySegment[1]) + 1
        for Gene, Sequence in Sequences.items():
            Identities[Gene] = float(len([Residues for Residues in zip(Sequences[ReferenceGene], Sequence)[UpStream:DownStream] \
                                              if Residues[0] == Residues[1]]))/len(Sequences[ReferenceGene][UpStream:DownStream])
    else:
        for Gene, Sequence in Sequences.items():
            Identities[Gene] = float(len([Residues for Residues in zip(Sequences[ReferenceGene], Sequence) \
                                              if Residues[0] == Residues[1]]))/len(Sequences[ReferenceGene])
    
    return Identities

class MutationObject:
    '''
    Class to handle various attributes of the mutation object.
    Assumes X#Z format. 
    '''

    def __init__(self, mutation):

        self.mutation = mutation 
        self.wild_type = mutation[:1]
        self.mutant   = mutation[-1:]
        self.residue  = int(mutation[1:-1])
        self.substitution = self.wild_type, self.mutant

def get_mutations(Arguments): #Classification, PathFile):
    '''
    Get the mutations from specified in the arguments file and pointed to in the paths file
    '''

    if Arguments.PostProcess == "Predict":
        return [Mutation.strip() for Mutation in file(get_file(Arguments.Mutations, Arguments.PathFile))]
    else:
        return dict([(Mutation.split()[0], float(Mutation.split()[1])) \
                         for Mutation in file(get_file(Arguments.Mutations, Arguments.PathFile))])
    
def get_annotation(Arguments):
    '''
    '''

    return dict([(int(Annotation.split()[0]), [Annotation.split()[1], float(Annotation.split()[2])]) \
                     for Annotation in file(get_file("Annotation", Arguments.PathFile))])
