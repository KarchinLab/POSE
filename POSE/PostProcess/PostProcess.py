'''
Get python modules
'''


'''
Get third-party modules
'''


'''
Get POSE modules
'''

from Validate import validate
from Predict import predict


def post_process(Arguments):
    ''' 
    This is simply a distibutor for all functions involved in 
    postprocessing POSE results from previous calculations, which 
    is required to make sense of most POSE output
    '''

    if Arguments.PostProcess == "Validate":
        validate(Arguments)

    if Arguments.PostProcess == "Predict":
        predict(Arguments)

    return 
