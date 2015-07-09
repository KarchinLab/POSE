'''
Get python modules
'''
from math import ceil
from random import shuffle
import cPickle

'''
Get third-party modules
'''

'''
Get POSE modules
'''
from MakePOSE import get_poses, pose_evaluator, epose_evaluator, make_pose_input, get_pose_scores
from GetData import get_mutations


def cross_validation(HoldoutMutations, Arguments):
   '''
   Take the data-splits dictionary from leave_some_out, and actually make the (e)POSEs for each 
   data split. 
   '''

   #Clustermode support if called
   Node = int(Arguments.MultiProcessMode[0]) 
   TotalNodes = int(Arguments.MultiProcessMode[1])

   Mutations, Sequences, ReferenceGene, Identities, ResidueBurial, Annotation = make_pose_input(Arguments)

   for CrossValidation, TestingMutations in HoldoutMutations.items()[Node::TotalNodes]:
      #Make a new dictionary without any of the testing mutations
      TrainingMutations = dict([(Mutation, Phenotype) for Mutation, Phenotype in Mutations.items() \
                                   if Mutation not in TestingMutations])
      
      POSEs = []
      Results = {}
      
      if Arguments.Mode == "POSE":
         for Trial in range(Arguments.Trials):
            POSEs.extend(get_poses(Arguments, TrainingMutations, Sequences, ReferenceGene, Identities, ResidueBurial, Annotation, \
                                      get_pose_scores, pose_evaluator))

      if Arguments.Mode == "ePOSE":
         for Trial in range(Arguments.Trials):
            POSEs.extend(get_poses(Arguments, TrainingMutations, Sequences, ReferenceGene, Identities, ResidueBurial, Annotation, \
                                      get_pose_scores, epose_evaluator))

      Results[tuple(TestingMutations)] = POSEs

      cPickle.dump(Results, open(Arguments.Filename + ".CrossValidation." + str(CrossValidation), "wb"), -1)

   return 


def leave_some_out(Arguments):
   '''
   Divide up the mutations before sending off to cross validation. 
   '''

   Mutations = get_mutations(Arguments)

   CrossValidations = int(ceil(len(Mutations.keys())/float(Arguments.LeaveSomeOut)))

   #For a standard POSE, the trick is keeping an as-balanced-as-possible spliting of the classes for each data split
   if Arguments.Mode == "POSE":
      Positive = [Mutation for Mutation, Phenotype in Mutations.items() if Phenotype]
      Negative = [Mutation for Mutation, Phenotype in Mutations.items() if not Phenotype]
   
      shuffle(Positive) #Get rid of bias that MIGHT be inherent in the original order of the mutations
      #Split as evenly as possible among the cross-validations
      Positive = dict([(CrossValidation, Positive[Mutation:len(Positive):CrossValidations]) \
                          for CrossValidation, Mutation in enumerate(range(CrossValidations))])
    
      shuffle(Negative) #Get rid of bias that MIGHT be inherent in the original order of the mutations
      #Split as evenly as possible among the cross-validations
      Negative = dict([(CrossValidation, Negative[Mutation:len(Negative):CrossValidations]) \
                          for CrossValidation, Mutation in enumerate(range(CrossValidations))])

      #only for the leave ONE out case do we do cases and controls in series
      if Arguments.LeaveSomeOut == 1:
         if len(Positive) < len(Negative):
            Postive = dict(zip(Positive.keys(), list(reversed(Positive.values()))))
         else:
            Negative = dict(zip(Negative.keys(), list(reversed(Negative.values()))))

      HoldoutMutations = dict([(CrossValidation, Positive[CrossValidation] + Negative[CrossValidation]) \
                                   for CrossValidation in Positive.keys()])

   #For ePOSEs it is simpler. Everyone is endophenotype positive. We therefore only have one class. 
   if Arguments.Mode == "ePOSE":
      Positive = Mutations.keys() 
      shuffle(Positive)  #Get rid of bias that MIGHT be inherent in the original order of the mutations
      HoldoutMutations = dict([(CrossValidation, Positive[Mutation:len(Positive):CrossValidations]) \
                                   for CrossValidation, Mutation in enumerate(range(CrossValidations))])
      
   cross_validation(HoldoutMutations, Arguments)

   return
   
   
