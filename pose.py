'''
Get python modules
'''
import pdb
import random
import cProfile 

'''
Get third-party modules
'''
from numpy import random as numdom

'''
Get POSE modules
'''
from POSE.MyPOSE import my_pose
from POSE.LeaveSomeOut import leave_some_out
from POSE.GetArguments import get_arguments
from POSE.MakePOSE import make_pose
from POSE.PostProcess.PostProcess import post_process

def main(Arguments):
    '''
    As of now, this routine exists just to call the appropriate 
    functions to run the desired MOCA mode.
    '''

    #Development/testing options
    if Arguments.Debug:
        pdb.set_trace()
    if Arguments.Seed:
        random.seed( Arguments.Seed )
        numdom.seed( Arguments.Seed ) #we need to make control numpy's randomness as well!
    
    if type(Arguments.Mutations) == list and len(Arguments.Mutations) == 1:
        Arguments.Mutations = Arguments.Mutations[0]

    if Arguments.Filename == "Default":
        Arguments.Filename =  Arguments.Mutations + "." + Arguments.Mode

    #What should I do with my POSE?
    POSEfnc = 0
    if Arguments.MakePOSE: #Make one
        make_pose(Arguments)
        POSEfnc += 1
    if Arguments.LeaveSomeOut:
        leave_some_out(Arguments)
        POSEfnc += 1
    if Arguments.PostProcess:
        post_process(Arguments)
        POSEfnc += 1
    if Arguments.MyPOSE: #.POSE == "MyPOSE": #Call a custom function that you added to MyPOSE.py
        my_pose(Arguments)
        POSEfnc += 1

    if not POSEfnc: 
        print "You ran POSE without specifying any POSE modes!!!"
        print "Exiting....."

if get_arguments().Profile:
    cProfile.run( "main(get_arguments())", sort="cumulative" )
else:
    main(get_arguments())
    
