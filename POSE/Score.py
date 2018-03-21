'''
Get python modules
'''
from collections import Counter
import itertools

'''
Get third-party modules
'''
from numpy import std, mean

'''
Get POSE modules
'''
import AminoAcidProperties as AA
from Statistics import cumulative_distribution_function as CDF

def molecular_weight_score(AminoAcids, AminoAcid):
    '''
    Make an identity-weighted histogram of molecular weights. Next, assume the histogram is normal. Take the molecular
    weight of your amino acid of interest and calculate the probability of observing an amino acid of that molecular
    weight in your normalized histogram.
    '''
    
    MolecularWeights = []
    for Key in AminoAcids:
        #We can't put fractions of an entry into our histogram bins; so, multiply by 1000 and round
        MolecularWeights.extend([AA.MolecularWeight[Key] for MolecularWeight in \
                                     range(int(round(1000*AminoAcids[Key])))])

    #print AminoAcids
    #print mean(MolecularWeights), std(MolecularWeights)
    #print AminoAcid, AA.MolecularWeight[AminoAcid]
    #print CDF(AA.MolecularWeight[AminoAcid] +5.0, mean(MolecularWeights), std(MolecularWeights)) - \
    #    CDF(AA.MolecularWeight[AminoAcid] - 5.0, mean(MolecularWeights), std(MolecularWeights))
    #Consider a molecular-weight window of 10 (i.e., 5 on either side of the weight of your aa of interest)
    return CDF(AA.MolecularWeight[AminoAcid] + 5.0, mean(MolecularWeights), std(MolecularWeights)) - \
        CDF(AA.MolecularWeight[AminoAcid] - 5.0, mean(MolecularWeights), std(MolecularWeights))

def hydropathy_index_score(AminoAcids, AminoAcid):
    '''
    Make an identity-weighted histogram of molecular weights. Next, assume the histogram is normal. Take the molecular
    weight of your amino acid of interest and calculate the probability of observing an amino acid of that molecular
    weight in your normalized histogram.
    '''
    
    HydropathyIndices = []
    for Key in AminoAcids:
        #We can't put fractions of an entry into our histogram bins; so, multiply by 1000 and round
        HydropathyIndices.extend([AA.HydropathyIndex[Key] for HydropathyIndex in \
                                     range(int(round(1000*AminoAcids[Key])))])

    #print AminoAcids
    #print mean(MolecularWeights), std(MolecularWeights)
    #print AminoAcid
    #print AA.MolecularWeight[AminoAcid]
    #print CDF(AA.MolecularWeight[AminoAcid], mean(MolecularWeights), std(MolecularWeights)) #- \
        #CDF(AA.MolecularWeight[AminoAcid] - 5.0, mean(MolecularWeights), std(MolecularWeights))
    #Consider a molecular-weight window of 10 (i.e., 5 on either side of the weight of your aa of interest)
    return CDF(AA.HydropathyIndex[AminoAcid] + 0.05, mean(HydropathyIndices), std(HydropathyIndices)) - \
        CDF(AA.HydropathyIndex[AminoAcid] - 0.05, mean(HydropathyIndices), std(HydropathyIndices))

def isoelectric_point_score(AminoAcids, AminoAcid):
    '''
    Make an identity-weighted histogram of molecular weights. Next, assume the histogram is normal. Take the molecular
    weight of your amino acid of interest and calculate the probability of observing an amino acid of that molecular
    weight in your normalized histogram.
    '''
    
    IsoelectricPoints = []
    for Key in AminoAcids:
        #We can't put fractions of an entry into our histogram bins; so, multiply by 1000 and round
        IsoelectricPoints.extend([AA.IsoelectricPoint[Key] for IsoelectricPoint in \
                                     range(int(round(1000*AminoAcids[Key])))])

    #print AminoAcids
    #print mean(MolecularWeights), std(MolecularWeights)
    #print AminoAcid
    #print AA.MolecularWeight[AminoAcid]
    #print CDF(AA.MolecularWeight[AminoAcid], mean(MolecularWeights), std(MolecularWeights)) #- \
        #CDF(AA.MolecularWeight[AminoAcid] - 5.0, mean(MolecularWeights), std(MolecularWeights))
    #Consider a molecular-weight window of 10 (i.e., 5 on either side of the weight of your aa of interest)
    return CDF(AA.IsoelectricPoint[AminoAcid] + 0.5, mean(IsoelectricPoints), std(IsoelectricPoints)) - \
        CDF(AA.IsoelectricPoint[AminoAcid] - 0.5, mean(IsoelectricPoints), std(IsoelectricPoints))

def amino_acid_score(AminoAcids, AminoAcid):
    '''
    Takes a dictionary of identity-weighted, amino-acid frequencies and returns that value for a requested 
    amino acid. Useful if dictionary represents a column in an MSA and the requested amino acid is a mutant
    or wild-type amino acid in a species of interest, from the same column. For instance, if a column
    has a well conserved alanine, and many species in the alignment have high identity relative to the species
    of interest, then an alanine will have a high score, and all other amino acids will have a low score. 
    '''

    return AminoAcids[AminoAcid]

def property_score(Properties, AminoAcid):
    '''
    Takes a dictionary of identity-weighted, amino-acid property frequencies and returns the sum of those 
    values for a paritular amino acid, normalized by the number of properties for that amino acid. Useful 
    if dictionary represents a column in an MSA and the requested amino acid is a mutant
    or wild-type amino acid in a species of interest, from the same column. For instance, if a column
    has a well conserved positive charge and donor capacity, and many species in the alignment have high 
    identity relative to the species of interest, then tyrosine will have a medium score becuase it can 
    donate a hydrogen, but is not positively charged. 
    '''

    return sum([Properties[Property] for Property in AA.Chemistry[AminoAcid]])/len(AA.Chemistry[AminoAcid])

def pose_score(MSA, mutation, Genes, Identities, ResidueBurial, Annotation, Arguments, GetProperties=False):
    '''
    Very basic score function that accounts for amino acid and biopyhsical properities conservation at 
    any requested column in an MSA. Provided an amino-acid substitution, it compares the wild type and mutant
    and returns a score based on the aforementioned conservations. 
    '''
    
    #You can't call "Structure" and score mutations not included in the structure!!! If so, return None for that residue
    if Arguments.Structure and mutation.residue not in ResidueBurial.keys(): return None
    
    #MolecularWeights = []
    NormalizationConstant = [] 
    WildType, Residue, Mutant = mutation.wild_type, mutation.residue - 1, mutation.mutant
    
    #Dictionaries with placeholders for each amino property or type
    Properties = dict([(Value, 0.0) for Values in AA.Chemistry.values() for Value in Values])
    AminoAcids = dict([(AminoAcid, 0.0) for AminoAcid in AA.Chemistry.keys()])
    
    #Build all the dictionaries/lists for the identity-weighted scores at the requested residue
    for Gene in Genes: 
        if MSA[Residue][Genes.index(Gene)] in AminoAcids:
            NormalizationConstant.append(Identities[Gene])
            AminoAcids[MSA[Residue][Genes.index(Gene)]] += Identities[Gene]
            for Property in AA.Chemistry[MSA[Residue][Genes.index(Gene)]]:
                Properties[Property] += Identities[Gene]

    #This tells you the fraction of sequences that did NOT have a X or - at the requested
    #position. X and - positions cannot be assessed and lead to decreasing confidence
    Confidence = float(len(NormalizationConstant))/len(Genes)
    
    if not Confidence or sum(NormalizationConstant) == 0:
        return None
    else:
        for AminoAcid in AminoAcids:
            AminoAcids[AminoAcid] = AminoAcids[AminoAcid]/sum(NormalizationConstant)
        
        for Property in Properties:
            Properties[Property] = Properties[Property]/sum(NormalizationConstant)
    '''
    WildTypeScore = amino_acid_score(AminoAcids, WildType) \
        + property_score(Properties, WildType) \
        + molecular_weight_score(AminoAcids, WildType)
    
    MutantScore = amino_acid_score(AminoAcids, Mutant) \
        + property_score(Properties, Mutant) \
        + molecular_weight_score(AminoAcids, Mutant)
    '''
    
    WildTypeScore = isoelectric_point_score(AminoAcids, WildType) \
        + hydropathy_index_score(AminoAcids, WildType) \
        + molecular_weight_score(AminoAcids, WildType) \
        + amino_acid_score(AminoAcids, Mutant)
    
    MutantScore = isoelectric_point_score(AminoAcids, Mutant) \
        + hydropathy_index_score(AminoAcids, Mutant) \
        + molecular_weight_score(AminoAcids, Mutant) \
        + amino_acid_score(AminoAcids, Mutant)
    
    if Arguments.Annotation:
        if mutation.residue in Annotation.keys():
            if Mutant in Annotation[mutation.residue][0]:
                WildTypeScore = WildTypeScore + Annotation[mutation.residue][1]
    
    if GetProperties:
        return molecular_weight_score(AminoAcids, WildType), molecular_weight_score(AminoAcids, Mutant), \
            hydropathy_index_score(AminoAcids, WildType), hydropathy_index_score(AminoAcids, Mutant), \
            isoelectric_point_score(AminoAcids, WildType), isoelectric_point_score(AminoAcids, Mutant)
   
    if Arguments.Structure:
        return (WildTypeScore - MutantScore)*ResidueBurial[mutation.residue]
    else:
        return WildTypeScore - MutantScore
