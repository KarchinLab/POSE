'''
Get python modules
'''

'''
Get third-party modulesz
'''

'''
Get POSE modules
'''

##########################################################
#
# Just an example to illustrate the MyPOSE function.
# 
# This exists so that users can combine pose functions in 
# different ways, and extend their usage here
#
# Under the "my_pose" function, put in as many new functions 
# as you want to call, and add them to this file. 
#
# See the functions below as examples
#
# To call, set --mypose fnc (from the command line) or 
# MyPOSE = fnc (from the Arguments file), where fnc is either
# "HelloWord" or "Goodbye" in the examples below
#
###########################################################

def hello_world(Arguments):
    '''
    Explain function here
    '''
   
    print "hello world"
   
    return 

def good_bye(Arguments):
    '''
    Explain function here
    '''
   
    print "good bye"
   
    return 

def my_pose(Arguments):
    
    if Arguments.MyPOSE == "HelloWorld":
        hello_world(Arguments)

    if Arguments.MyPOSE == "GoodBye":
        good_bye(Arguments)

    return
