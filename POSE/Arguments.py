'''
Get python modules
'''
from argparse import ArgumentParser as argument_parser
from sys import argv
from ast import literal_eval
import re
import os

'''
Get third-party modules
'''

'''
Get POSE modules
'''

def get_argument_file(Arguments):
    '''
    This function is used to read in an agrument file, instead of reading arguments
    from the command line. First, all defaults are set from the get_arguments function, 
    then overide defaults with any command that appears in the Arguments file. 
    This mode will run if there are either no commands, in which case looks for ./Arguments--
    if it finds this file, it will read it in. Alternatively, you can tell MOCA the path/name
    of the arguments file from the command line with the --argumentfile (-a) command.
    '''

    class ArgumentObject:

        def __init__(self, Arguments):

            for Argument in Arguments:
                try: #first see if the argument is number
                    Value = float(Argument[1])
                except ValueError:
                    Value = Argument[1]
                if Argument[1].isdigit(): #specifically, is it an int?
                    vars(self)[Argument[0]] = int(Argument[1])
                elif Argument[1] == "True" or Argument[1] == "False": #should we convert to bool?
                    vars(self)[Argument[0]] = literal_eval(Value)
                elif type(Value) == float:
                    vars(self)[Argument[0]] = Value
                else: #do we need to convert to a python list?
                    if len(re.sub(r"[?]", " ",  Value).split()) > 1:
                        vars(self)[Argument[0]] = re.sub(r"[?]", " ",  Value).split()
                    else:
                        vars(self)[Argument[0]] = Value

    #Make the argument-value pairs
    UserArguments = [(Argument.split()[0], " ".join(Argument.split()[2:]).replace(",", "")) \
                     for Argument in file(Arguments.ArgumentFile) \
                     if Argument.strip() and Argument.split()[1] == "="]
    
    #Now, overide any default arguments with those specified in the arguments file
    for Key, Value in ArgumentObject(UserArguments).__dict__.items():
        Arguments.__dict__[Key] = Value
    
    return Arguments

def get_arguments():
    '''
    Read in the POSE defaults. Override defaults with from the command line or arguments file
    See UsersManual.pdf for more detail on commands
    '''

    parser = argument_parser()
    parser.add_argument("-a", "--argumentfile", dest="ArgumentFile", default="Arguments", \
                            help="Use an argument file rather than the command line. Supply the name of the argument file. Note: POSE automatically looks for './Arguments'.")
    parser.add_argument("-b", "--seed", dest="Seed", default=0, type=int, \
                            help="Supply a 'seed' (and supply a long int) when output reproducility is required (controls stochastic components of algorithm)")
    parser.add_argument("-c", "--correlated", dest="Correlated", action="store_true", default=True, \
                            help="Do increasing POSE score correlate with the value of the measurement, or are these numbers anticorrelated?")
    parser.add_argument("-d", "--domain", dest="Domain", nargs=2, default=[0, 0], \
                            help="User-defined range of amino acids to be considered at either the initial Blastp and alignment, and during (e)POSE derivation. Example: 'Domain = 305, 678'. The real reason to use this argument is to isolate distinct evolutionary modules within you protein of interest.")
    parser.add_argument("-e", "--email", dest="Email", default=None, \
                            help="Give me the name you used to point to mutation files in the paths file. See the Paths and Arguments files from the tutorial for example POSE and ePOSE usage")
    parser.add_argument("-f", "--filename", dest="Filename", default="Default", \
                            help="Input/output filname prefix for reading/writing files. By default, filename prefix is the string provided to the mutation command + '.' Mode (e.g., CFTR.POSE")
    parser.add_argument("-g", "--geneidentifier", dest="GeneIdentifier", action="store_true", default=None, \
                            help="Perhaps more appropriately a 'protein identifier', this can be either a hugo gene symbol (e.g., PIK3CA) or a uniprot identifier (e.g., P42336)")
    #h = help
    parser.add_argument("-i", "--optimizationparameters", dest="OptimizationParameters", nargs=3, default=[100, 0.01, 10], \
                            help="The first param specifies the frequency to repopulate the sequence pool with top predicting POSEs, the second param specificies the percent of top POSEs to repopulate the pool, and the third is the number to write to disk for each trial. Defaults specify that every 100 optimizations the top-performing 1%% of POSEs go back into the sequence pool, and when a whole trial is done, write the top 10 POSEs to disk.")
    parser.add_argument("-j", "--preprocess", dest="PreProcess", action="store_true", default=None, \
                            help="FILL IN")
    parser.add_argument("-k", "--mutations", dest="Mutations", nargs="*", \
                            help="Give me the name you used to point to mutation files in the paths file. See the Paths and Arguments files from the tutorial for example POSE and ePOSE usage")
    #parser.add_argument("-k", "--identitysegment", dest="IdentitySegment", nargs=2, default=[0, 0], \
     #                       help="Sometimes you don't want the entire protein sequence to be included when determining the sequence weights. For instance, you might be scoring mutations from a evolutionarily modular domain, so some of the full sequences in the initial pool might be irrelevant outside of this domain. Accounts for python slicing (i.e., enter the exact resiues you want 'upstream', 'downstream')")
    parser.add_argument("-l", "--leavesomeout", dest="LeaveSomeOut", default=0, type=int, \
                            help="Specificy the number to leave out for a leave-some-out cross-validation")
    parser.add_argument("-m", "--mode", dest="Mode", default="PreProcess", 
                        help="'PreProcess' if you want POSE to prepare the initial sequence pool (i.e., input MSA) for you. '(e)POSE' to make either an ePOSE or POSE. 'PostProcess' for analyzing results. See POSE user's manual for all options (default=PreProcess)")
    parser.add_argument("-n", "--maxrandomsequences", dest="MaxRandomSequences", default=150, type=int, \
                            help="Maxium number of random sequences sampled during a single iteration of POSE construction")
    parser.add_argument("-o", "--optimizations", dest="Optimizations", default=1000, type=int, \
                            help="How many many random sequences samples to rank for a single trial (t)")
    parser.add_argument("-p", "--pathfile", dest="PathFile", default="Paths", \
                            help="Show me where the path file is. Default is './Paths'")
    parser.add_argument("-q", "--chain", dest="Chain", default=" ", \
                            help="If structure is called, should I be looking at a particular chain in the PDB file. Note: Burial calculations consider all atoms in the PDB, but I need to know which residues to consider for scoring purposes")
    parser.add_argument("-r", "--makepose", dest="MakePOSE", action="store_true", default=False, \
                            help="Default mode for constructing phenotype-optimized sequences ensembles")
    parser.add_argument("-s", "--structure", dest="Structure", action="store_true", default=False, \
                            help="Provide PDB-formatted structure file to calculate and consider residue burial for scoreing mutations")
    parser.add_argument("-t", "--trials", dest="Trials", default=1, type=int, \
                            help="How many times should I start POSE optimization from the beginning. This is important because each start can lead to different local minima")
    parser.add_argument("-u", "--multiprocessmode", dest="MultiProcessMode", nargs=2, default=[0, 1], \
                            help="If distributing jobs accross multiple processors, the first option determines which node is running the job, and the second argument determines how many nodes are running")
    parser.add_argument("-v", "--annotation", dest="Annotation", action="store_true", default=False, \
                            help="Provid an annotation file if desired. e.g.: the line '231 HGF 0.5' would mean if residue 231 is mutated to either an H, G, or F, then further penalize the mutation by 0.5 score units. You can have as many lines as you like in your annotation file.")
    parser.add_argument("-w", "--postprocess", dest="PostProcess", action="store_true", default=False, \
                            help="Choose a method for postproccessing your POSE results")

    #Post-processing options get capital letters

    #Developer options
    parser.add_argument("--mypose", dest="MyPOSE", action="store_true", default=False, 
                        help="Build your own MOCA module! See MyMOCA.py")
    parser.add_argument("--profile", dest="Profile", action="store_true", default=False,
                        help="Profile to find inefficient code blocks")
    parser.add_argument("--debug", dest="Debug", action="store_true", default=False,
                        help="run in python debug mode 'pdb'")
    
    Arguments = parser.parse_args()
    
    #Overwrite  any defaults from the Arguments file
    if os.path.isfile(Arguments.ArgumentFile): Arguments = get_argument_file(Arguments)
    
    return Arguments  
