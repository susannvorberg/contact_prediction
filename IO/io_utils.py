import os
import json

def read_json_from_mat(matfile):
    '''
        Read the specified keys from the json data
        line with json data must start with #META
    :param matfile: contact matrix file
    :return: return dict of meta data
    '''

    if not os.path.exists(matfile):
        raise FileNotFoundError("Specified matfile does not exist: " + str(matfile))

    meta={}

    with open(matfile, 'r') as f:
        for line in f:
            if '#META' in line:
                meta = json.loads(line)
            else:
                print(str(matfile) + "does not contain META info. (Line must start with #META!)")

    return meta