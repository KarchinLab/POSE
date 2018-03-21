'''
Get python modules
'''
import cPickle
from collections import Counter
'''
Get third-party modules
'''
from numpy import mean, std

'''
Get POSE modules
'''
from POSE.MakePOSE import make_pose_input, get_pose_scores
from POSE.SystemCommands import ls
from POSE.Statistics import ROC, linear_regression, correlation, correlation_pvalue

def pose_prediction(Arguments):
    '''
    Score some new mutations with the POSE(s) you already created. First we need to determine the optimal POSE
    score cutoff for predicting a mutants phenotype; this is done using the mutations used to create the POSE in
    the first place. Next, we score the mutations and apply said cutoffs to make predictions. 
    '''

    #First get the cutoff from the mutations originally used to make the POSE
    Mutations, Sequences, ReferenceGene, Identities, ResidueBurial, Annotation = make_pose_input(Arguments)
    #These mutations are the ones we actually need.

    Mutations, POSEs = zip(*[(cPickle.load(open(File, "rb"))[0], cPickle.load(open(File, "rb"))[1][0]) \
                             for File in ls("./") if Arguments.Filename in File])
    Mutations, Phenotypes = Mutations[0].keys(), Mutations[0].values()

    Scores = []
    for POSE in POSEs:
        MSA = zip(*[Sequences[Gene] for Gene in POSE]) #Turn the list of genes into the MSA (sequence ensemble)
        Scores.append(get_pose_scores(MSA, list(POSE), Identities, Mutations, ResidueBurial, Annotation, Arguments))

    Scores = [mean(Score) for Score in zip(*Scores)]

    AUC, Metrix = ROC(Scores, Phenotypes, Arguments.Correlated)
    Optimal = sorted(Metrix, key=lambda Metric: Metric.Sensitivity + Metric.Specificity, reverse=True)[0]
    
    #Now we are ready to get the new mutations. 
    Mutations, Sequences, ReferenceGene, Identities, ResidueBurial, Annotation = make_pose_input(Arguments)

    Scores = []
    for POSE in POSEs:
        MSA = zip(*[Sequences[Gene] for Gene in POSE]) #Turn the list of genes into the MSA (sequence ensemble)
        Scores.append(get_pose_scores(MSA, list(POSE), Identities, Mutations, ResidueBurial, Annotation, Arguments))

    Scores = dict(zip(Mutations, zip(*Scores)))

    print "Mutation, Mean POSE score, Standard deviation, Predicted phenotype"
    for Mutation in Mutations:
        if Arguments.Correlated and mean(Scores[Mutation]) >= Optimal.CutOff:
            print Mutation, "\t", "\t", round(mean(Scores[Mutation]), 3), "\t", "\t", round(std(Scores[Mutation]), 3), "\t", "\t", "Positive"

        if Arguments.Correlated and mean(Scores[Mutation]) < Optimal.CutOff:
            print Mutation, "\t", "\t", round(mean(Scores[Mutation]), 3), "\t", "\t", round(std(Scores[Mutation]), 3), "\t", "\t", "Negative"

        if not Arguments.Correlated and mean(Scores[Mutation]) <= Optimal.CutOff:
            print Mutation, "\t", "\t", round(mean(Scores[Mutation]), 3), "\t", "\t", round(std(Scores[Mutation]), 3), "\t", "\t", "Positive"

        if not Arguments.Correlated and mean(Scores[Mutation]) > Optimal.CutOff:
            print Mutation, "\t", "\t", round(mean(Scores[Mutation]), 3), "\t", "\t", round(std(Scores[Mutation]), 3), "\t", "\t", "Negative"

    print " "
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    print "Cutoff used for separating the classes (determined when POSEs were initially created!!!) =", round(Optimal.CutOff, 2)
    print " "
    print "Performance acheived on the mutations used to train the initial POSE:"
    print " "
    print "\t", "Area under the ROC curve =", round(AUC, 2)
    print "\t", "Accuracy =", Optimal.Accuracy
    print "\t", "Sensitivity =", Optimal.Sensitivity
    print "\t", "Specificity =", Optimal.Specificity
    print "\t", "Negative predictive value =", Optimal.NPV
    print "\t", "Positive predictive value =", Optimal.PPV
    print "\t", "Mathews correlation coefficient =", Optimal.MCC
    
    return

def epose_prediction(Arguments):
    '''
    Score some new mutations with the ePOSE(s) you already created. First we'll take our ePOSEs and the mutations originally used 
    create them to derive the linear equation necessary to convert POSE scores to predicted endophenotypes. Finally, we score the new
    mutations and use said linear equation to supply predictions. 
    '''

    #First get the cutoff from the mutations originally used to make the ePOSE
    Mutations, Sequences, ReferenceGene, Identities, ResidueBurial, Annotation = make_pose_input(Arguments)
    #These mutations are the ones we actually need.
    Mutations, ePOSEs = zip(*[(cPickle.load(open(File, "rb"))[0], cPickle.load(open(File, "rb"))[1][0]) \
                             for File in ls("./") if Arguments.Filename in File])
    Mutations, Endophenotypes = Mutations[0].keys(), Mutations[0].values()

    LinearParameters = []
    ePOSEScores = []
    for ePOSE in ePOSEs:
        MSA = zip(*[Sequences[Gene] for Gene in ePOSE])  #Turn the list of genes into the MSA (sequence ensemble)
        
        Scores = get_pose_scores(MSA, list(ePOSE), Identities, Mutations, ResidueBurial, Annotation, Arguments)
        ePOSEScores.append(Scores)
        a, b, RR = linear_regression(Scores, Endophenotypes)
        
        LinearParameters.append((a, b))

    ePOSEScores = [mean(Score) for Score in zip(*ePOSEScores)]

    #Now we are ready to get the new mutations. 
    Mutations, Sequences, ReferenceGene, Identities, ResidueBurial, Annotation = make_pose_input(Arguments)
    Scores = []
    
    for ePOSE in ePOSEs:
        MSA = zip(*[Sequences[Gene] for Gene in ePOSE]) #Turn the list of genes into the MSA (sequence ensemble)
        Scores.append(get_pose_scores(MSA, list(ePOSE), Identities, Mutations, ResidueBurial, Annotation, Arguments))
    
    Scores = dict(zip(Mutations, zip(*Scores)))
    
    print "Mutation, Mean POSE score, Standard deviation, Predicted Endophenotype, Standard deviation"
    for Mutation in Mutations:
        PredictedEndophenotypes = [a*mean(Scores[Mutation]) + b for a, b in LinearParameters] 
        print Mutation, "\t", "\t", round(mean(Scores[Mutation]), 3), "\t", "\t", round(std(Scores[Mutation]), 3), \
            "\t", "\t", "\t", round(mean(PredictedEndophenotypes), 3), "\t", "\t", "\t", round(std(PredictedEndophenotypes), 3)

    print " "
    print "~~~~~~~~~~ Performance Metrics ~~~~~~~~~~~~~"
    print " "
    print "R-squared =", linear_regression(ePOSEScores, Endophenotypes)[2]
    print "Pearson correlation =", correlation(ePOSEScores, Endophenotypes)
    print "P-value =", correlation_pvalue(ePOSEScores, Endophenotypes)
        
    return 

def predict(Arguments):
    '''
    With your (e)POSEs, created using MakePOSE, in the current-working directory, point POSE to your list
    of new mutations to score. 
    '''
    
    if Arguments.Mode == "POSE": pose_prediction(Arguments)

    if Arguments.Mode == "ePOSE": epose_prediction(Arguments)
    
    return

