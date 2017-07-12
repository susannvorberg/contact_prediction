
import os
import numpy as np
import json
import pandas as pd
from scipy.stats import norm, multivariate_normal



def gaussian_mixture_density_2d(x, y, weights, means_ab, means_cd, covMats, log=False):
    if (log):
        sum_density = np.log(1e-323)
    else:
        sum_density = 0
    for component in range(len(weights)):
        if (log):
            sum_density = logAdd(sum_density, np.log(weights[component]) + multivariate_normal.logpdf([x, y], [
                means_ab[component], means_cd[component]], covMats[component]))
        else:
            sum_density += weights[component] * multivariate_normal.pdf([x, y],
                                                                        [means_ab[component], means_cd[component]],
                                                                        covMats[component])
    return sum_density

def gaussian_mixture_density(x, weights, means, sd):

        sum_density = 0
        for component in range(len(weights)):
           sum_density +=  weights[component] * norm.pdf(x, means[component], sd[component])
        return sum_density


def read_optimization_log_file(optimization_log_file):

    if not os.path.exists(optimization_log_file):
        print("Optimization log file {0} does not exist!".format(optimization_log_file))
        return  pd.DataFrame({})

    try:
        with open(optimization_log_file) as json_data:
            json_string = json.load(json_data)
    except Exception as e:
        print("Optimization log file {0} is not in json format!:{1}".format(optimization_log_file, e))
        return  pd.DataFrame({})

    log_df = pd.read_json(json_string, orient='split')
    return log_df

def read_parameter_file(parameter_file):

    parameters={}

    #read parameters
    if not os.path.exists(parameter_file):
        print("Parameter file {0} does not exist!".format(parameter_file))
        return parameters

    try:
        with open(parameter_file, 'r') as f:
            parameters = json.load(f)
    except IOError:
        print("Cannot open {0} in json format!".format(parameter_file))
        return parameters

    return parameters

def read_settings_file(settings_file):

    settings = {}

    #read settings
    if not os.path.exists(settings_file):
        print("Settings file {0} does not exist!".format(settings_file))

    with open(settings_file, 'r') as f:
        settings = json.load(f)


    return settings

