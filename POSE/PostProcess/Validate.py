'''
Get python modules
'''
import cPickle


'''
Get third-party modules
'''
from numpy import mean, std

'''
Get POSE modules
'''
from POSE.SystemCommands import ls
from POSE.MakePOSE import make_pose_input, get_pose_scores
from POSE.Statistics import ROC, linear_regression, correlation, correlation_pvalue

def validate_pose(Arguments):
    '''
    Validate the POSE. See the validate function (below) for detailed description. 
    '''

    CrossValidations = [cPickle.load(open(File, "rb")) for File in ls("./") if Arguments.Filename + ".CrossValidation." in File]

    Mutations, Sequences, ReferenceGene, Identities, ResidueBurial, Annotation = make_pose_input(Arguments)
    
    Predictions = {}
    for CrossValidation in CrossValidations:
        TestingMutations, POSEs = list(CrossValidation.keys()[0]), CrossValidation.values()[0]

        TestScores = []
        for POSE in POSEs:
            MSA = zip(*[Sequences[Gene] for Gene in POSE]) #Turn the list of genes into the MSA (sequence ensemble)
            TestScores.append(get_pose_scores(MSA, list(POSE), Identities, TestingMutations, ResidueBurial, Annotation, Arguments))
        
        #If you had more than one trial in you MakePOSE, you have more than one score associated with each mutation. Here we collect em
        for Mutation, Scores in dict(zip(TestingMutations, zip(*TestScores))).items():
            #Now we get the mean and standard deviation of scores associated with each mutation (relevant if trial > 1 in you initial MakePOSE).
            Predictions[Mutation] = [round(mean(Scores), 3), round(std(Scores), 3)]

            
    print ", ".join(["Mutation", "'Known' phenotype", "POSE score", "Std dev (if multiple 'trials' were used)"])
    for Mutation in Mutations.keys():
        print Mutation + "\t" + "\t" + ["Positive" if Mutations[Mutation] else "Negative"][0] \
            + "\t" + "\t ".join(map(str, Predictions[Mutation]))
        
    print " "
    print "~~~~~~~~~~ Performance Metrics ~~~~~~~~~~~~~"
    Scores, Phenotypes = zip(*[(Predictions[Mutation][0], int(Mutations[Mutation])) for Mutation in Mutations.keys()])
    AUC, Metrix = ROC(Scores, Phenotypes, Arguments.Correlated)
    Optimal = sorted(Metrix, key=lambda Metric: Metric.Sensitivity + Metric.Specificity, reverse=True)[0]
    
    print "Optimal cutoff for separating the classes =", round(Optimal.CutOff, 2)
    print " "
    print "Performance using optimal cutoff:"
    print " "
    print "\t", "Area under the ROC curve =", round(AUC, 2)
    print "\t", "Accuracy =", Optimal.Accuracy
    print "\t", "Sensitivity =", Optimal.Sensitivity
    print "\t", "Specificity =", Optimal.Specificity
    print "\t", "Negative predictive value =", Optimal.NPV
    print "\t", "Positive predictive value =", Optimal.PPV
    print "\t", "Mathews correlation coefficient =", Optimal.MCC
    
    return

def validate_epose(Arguments):
    '''
    Validate the ePOSE. See the validate function (below) for detailed description. 
    '''

    CrossValidations = [cPickle.load(open(File, "rb")) for File in ls("./") if Arguments.Filename + ".CrossValidation." in File]

    Mutations, Sequences, ReferenceGene, Identities, ResidueBurial, Annotation = make_pose_input(Arguments)

    Predictions = {}
    for CrossValidation in CrossValidations:
        TestingMutations, ePOSEs = list(CrossValidation.keys()[0]), CrossValidation.values()[0]

        TrainingMutations = dict([(Mutation, Endophenotype) for Mutation, Endophenotype in Mutations.items() \
                                      if Mutation not in TestingMutations])
        
        TestScores = []
        PredictedEndophenotypes = [] #this is where will put ePOSE scores that are converted to their corresponding endophenotype prediction
        for ePOSE in ePOSEs:
            MSA = zip(*[Sequences[Gene] for Gene in ePOSE])  #Turn the list of genes into the MSA (sequence ensemble)
            
            test_scores = get_pose_scores(MSA, list(ePOSE), Identities, TestingMutations, ResidueBurial, Annotation, Arguments)
            TestScores.append(test_scores)

            TrainScores = get_pose_scores(MSA, list(ePOSE), Identities, TrainingMutations.keys(), ResidueBurial, Annotation, Arguments)
            a, b, RR = linear_regression(TrainScores, TrainingMutations.values())
            PredictedEndophenotypes.append([a*x + b for x in test_scores])
            
        PredictedEndophenotypes = dict(zip(TestingMutations, zip(*PredictedEndophenotypes)))

        #If you had more than one trial in you MakePOSE, you have more than one score associated with each mutation. Here we collect em
        for Mutation, Scores in dict(zip(TestingMutations, zip(*TestScores))).items():
            #Now we get the mean and standard deviation of scores associated with each mutation (relevant if trial > 1 in you initial MakePOSE).
            Predictions[Mutation] = [round(mean(Scores), 3), round(std(Scores), 3), \
                                         round(mean(PredictedEndophenotypes[Mutation]), 3), round(std(PredictedEndophenotypes[Mutation]), 3)]
        
    Scores = [Predictions[Mutation][0] for Mutation in Mutations.keys()] #make sure scores are ordered just like Mutations.values()
    print ", ".join(["Mutation", "Measurement", "Score", "Std dev", "Prediction", "Std dev"])
    for Mutation in Mutations.keys():
        print Mutation + "\t" + "\t" + str(Mutations[Mutation]) + "\t" + "\t ".join(map(str, Predictions[Mutation]))

    print " "
    print "~~~~~~~~~~ Performance Metrics ~~~~~~~~~~~~~"
    print " "
    print "R-squared =", linear_regression(Scores, Mutations.values())[2]
    print "Pearson correlation =", correlation(Scores, Mutations.values())
    print "P-value =", correlation_pvalue(Scores, Mutations.values())

    return

def validate(Arguments):
    '''
    Get performance for your leave-some-out cross-validation (e)POSE run. The dictionary you load is holdout-mutation keys
    and the sequence ensembles selected for the corresponding training mutations (i.e., everything but the holdout mutations). Now
    we'll score the holdout mutations and see how we did. Because phenotypes are dichotomous (binary) and endophenotypes are 
    continuous, we need different score functions for each. 
    '''

    if Arguments.Mode == "POSE": validate_pose(Arguments)

    if Arguments.Mode == "ePOSE": validate_epose(Arguments)

    return
