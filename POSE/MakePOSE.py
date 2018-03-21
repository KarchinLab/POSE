'''
Get python modules
'''
from random import sample, randint
import cPickle

'''
Get third-party modules
'''

'''
Get POSE modules
'''
from GetData import get_sequences, get_reference_gene, get_identities, get_mutations, MutationObject, \
    get_annotation
from Score import pose_score
from Statistics import ROC, linear_regression
from Structure import normalized_residue_burial

def get_random_genes(GenePool, MaxRandomSequences):
    '''
    During each optimization iteration of POSE construction, a random group of genes is pulled from the 
    gene pool. This function gets those genes, which can be as few as 1 and a many as the max, which is 
    set in the command arguments/arguments file. 
    '''
    
    return sample(GenePool, min(len(GenePool), randint(1, MaxRandomSequences)))

def get_pose_scores(MSA, Genes, Identities, Mutations, ResidueBurial, Annotation, Arguments):
    '''
    This basically calls the default pose score function, for each mutation one at a time. It returns a list
    of all mutations for which the provided MSA facilitated scoring; a mutation cannot be scored if every gene, 
    at that column in the alignment, had an 'X' or a '-'. If you prefer a different score function, you can 
    provide one to the make_pose function. All you need to do is edit the MyPOSE script, and call MyPOSE rather 
    than MakePOSE.
    '''
    
    Scores = [pose_score(MSA, MutationObject(Mutation), Genes, Identities, ResidueBurial, Annotation, Arguments) \
                      for Mutation in Mutations]
    
    return filter(lambda Score: Score != None, Scores)

def pose_evaluator(SequenceEnsembles, TopPOSEs, Arguments):
    '''
    The evaluator assesses the performance/predictive value of each potential POSE. Input is a dictionary 
    of potential POSE keys, each with a list of scores and a list of phenotypes (ones and zeros)
    associated with each of the potential POSEs. The default evaluator uses the area under the ROC curve 
    (AUC). If you prefer a different evaluation metric, you can provide one to the make_pose function. All 
    you need to do is edit the MyPOSE script, and call MyPOSE rather than MakePOSE.
    '''

    Ensembles = SequenceEnsembles.keys()
    Scores, Phenotypes = zip(*SequenceEnsembles.values())
    Phenotypes = [map(int, Phenotype) for Phenotype in Phenotypes] #change floats to ints for the hell of it (0.0 = 0; 1.0 = 1)
    
    EvaluatedEnsembles = {}
    for Ensemble in Ensembles:
        AUC, Metrix = ROC(Scores[Ensembles.index(Ensemble)], Phenotypes[Ensembles.index(Ensemble)], Arguments.Correlated)
        try:
            Optimal = sorted(Metrix, key=lambda Metric: Metric.Sensitivity + Metric.Specificity, reverse=True)[0]
            EvaluatedEnsembles[Ensemble] = AUC #could do something like Optimal.Sensitivity + Optimal.Specificity
        except IndexError:
            pass

    POSEs = []
    for Ensemble, Performance in sorted(EvaluatedEnsembles.items(), key=lambda Items: Items[1], reverse=True)[:TopPOSEs]:
        POSEs.append(Ensemble)
    
    return POSEs

def epose_evaluator(SequenceEnsembles, TopPOSEs, Arguments):
    '''
    This evaluator assesses the performance of each potential ePOSE. Input is a dictionary 
    of potential ePOSE keys, each with a list of scores and corresponding endophenotypes
    associated with each of the potential ePOSEs. This default evaluator uses the R-squared 
    of scores and endophenotypes. If you prefer a different evaluation metric, you can provide 
    one to the make_pose function. All you need to do is edit the MyPOSE script, and call 
    MyPOSE rather than MakePOSE.
    '''
    
    Ensembles = SequenceEnsembles.keys()
    Scores, EndoPhenotypes = zip(*SequenceEnsembles.values())

    EvaluatedEnsembles = {}
    for Ensemble in Ensembles:
        a, b, RR = linear_regression(Scores[Ensembles.index(Ensemble)], EndoPhenotypes[Ensembles.index(Ensemble)])

        #The slope and correlation should match (i.e., we don't want positive slope for an anti-correlated interaction) 
        if a < 0.0 and not Arguments.Correlated:
            EvaluatedEnsembles[Ensemble] = [a, b, RR]

        #The slope and correlation should match (i.e., we don't want negative slope for an correlated interaction) 
        if a > 0.0 and Arguments.Correlated:
            EvaluatedEnsembles[Ensemble] = [a, b, RR]
    
    POSEs = []
    for Ensemble, Performance in sorted(EvaluatedEnsembles.items(), key=lambda Items: Items[1][2], reverse=True)[:TopPOSEs]:
        POSEs.append(Ensemble)
    
    return POSEs

def make_pose_input(Arguments):
    '''
    Gather all the data required as input for (e)POSE derivation. 
    '''

    Mutations = get_mutations(Arguments) #Dictionary of amino acid substitutions and binary pheontypes or endophenotypes
    Sequences = get_sequences(Arguments) #Load the fasta formatted sequence file
    ReferenceGene = get_reference_gene(Arguments) #Gene your scoring
    Identities = get_identities(ReferenceGene, Sequences, Arguments) #%ID of all sequences relative to ref

    #Inititialize burial, and populate if called
    ResidueBurial = {}
    if Arguments.Structure: ResidueBurial = normalized_residue_burial(Arguments)

    #Inititialize annotate, and populate if called
    Annotation = {}
    if Arguments.Annotation: Annotation = get_annotation(Arguments)

    return Mutations, Sequences, ReferenceGene, Identities, ResidueBurial, Annotation

def get_poses(Arguments, Mutations, Sequences, ReferenceGene, Identities, ResidueBurial, Annotation, \
                  score_function, evaluation_method):
    '''
    The main function for making phenotype-optimized sequence ensembles (POSEs), using the default score
    function and performance evaulation method. If you prefer a different score function and/or
    evaluation method, you can provide one to the make_pose function. All you need to do is edit the MyPOSE 
    script, and call MyPOSE rather than MakePOSE.
    '''

    '''
    How often to repopulate gene pool, what percent of top-performing POSEs to add to pool, 
    and what percent to save at the end
    '''
    RepopulateFrequency = int(Arguments.OptimizationParameters[0])
    POSEsToRepopulate = float(Arguments.OptimizationParameters[1])
    POSEsToWrite = int(Arguments.OptimizationParameters[2])

    
    #for Trial in range(Arguments.Trials)[Node::TotalNodes]: #Number of starts from initial pool. 
    GenePool = Sequences.keys()
    SequenceEnsembles = {}
    for Optimization in range(Arguments.Optimizations):

        Genes = get_random_genes(GenePool, Arguments.MaxRandomSequences)
        
        #Each tuple is a column in the alignment with residues for all rand genes
        MSA = zip(*[Sequences[Gene] for Gene in Genes])

        #Score the mutations. Mutation.keys() are the amino acid substitutions
        Scores = score_function(MSA, Genes, Identities, Mutations.keys(), ResidueBurial, Annotation, Arguments)
        
        #Make sure all mutations were scorable with the most recent sequence ensemble. Mutations.values() are the (endo)phenotypes
        if len(Scores) == len(Mutations.values()): SequenceEnsembles[tuple(Genes)] = [Scores, Mutations.values()]
        
        if not Optimization%RepopulateFrequency and Optimization > 0:
            POSEs = evaluation_method(SequenceEnsembles, int(POSEsToRepopulate*len(SequenceEnsembles.keys())), \
                                          Arguments)
            [GenePool.append(Gene) for POSE in POSEs for Gene in POSE]
            
    return evaluation_method(SequenceEnsembles, POSEsToWrite, Arguments)
            
def make_pose(Arguments):
    '''
    Make the POSE or endophenotyp POSE (ePOSE)
    '''

    Mutations, Sequences, ReferenceGene, Identities, ResidueBurial, Annotation = make_pose_input(Arguments)

    #Clustermode support if called
    Node = int(Arguments.MultiProcessMode[0]) 
    TotalNodes = int(Arguments.MultiProcessMode[1]) 

    for Trial in range(Arguments.Trials)[Node::TotalNodes]:
        if Arguments.Mode == "POSE":
            POSE = get_poses(Arguments, Mutations, Sequences, ReferenceGene, Identities, ResidueBurial, Annotation, get_pose_scores, pose_evaluator)
        if Arguments.Mode == "ePOSE":
            POSE = get_poses(Arguments, Mutations, Sequences, ReferenceGene, Identities, ResidueBurial, Annotation, get_pose_scores, epose_evaluator)

        cPickle.dump([Mutations, POSE], open(Arguments.Filename + "." + str(Trial), "wb"), -1)

    return
