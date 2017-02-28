import numpy as np
import json
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

def init_by_default(initial_weights, initial_means, initial_precision, sigma, fixed_zero_component=True):


    parameters = {}


    parameters['weight_contact_0'] = [initial_weights[0]]  ##* 400
    parameters['weight_bg_0'] = [initial_weights[0]]  ##* 400
    parameters['mu_0'] = [initial_means[0]] * 400
    parameters['prec_0'] = [initial_precision[0]] * 400

    start_range = 0
    if(fixed_zero_component):
        start_range =1

    for component in range(start_range, len(initial_weights)):

        parameters['weight_contact_' + str(component)] = [initial_weights[component]]  ##* 400
        parameters['weight_bg_' + str(component)] = [initial_weights[component]]  ##* 400
        parameters['mu_' + str(component)] = np.random.normal(loc=initial_means[component],
                                                              scale=0.05,
                                                              size=400).tolist()

        initial_diagonal = np.random.normal(
            loc=initial_precision[component],
            scale=initial_precision[component] / 10,
            size=400
        )

        precMat = np.zeros((400, 400))
        np.fill_diagonal(precMat,initial_diagonal)


        if (sigma == 'diagonal'):
            parameters['prec_' + str(component)] = np.diag(precMat).tolist()
        if (sigma == 'full'):
            parameters['prec_' + str(component)] = precMat[np.tril_indices(400)].tolist()

    return parameters


def transform_parameters(parameters_in, weights=True, mean=False, prec=True, back=False):
    """
    parameters: Pandas dataframe with column names
    weights:    will be transformed with softmax
    mean:       will be transformed with exp
    prec:       will be transformed with exp
    back:       will regain original parameters
    """

    transformed_parameters = parameters_in.copy()

    if (not back):

        # weights will be represented as softmax
        if (weights):
            weight_names_bg = [parameter_name for parameter_name in parameters_in.keys() if
                               "weight_bg" in parameter_name]
            weight_names_contact = [parameter_name for parameter_name in parameters_in.keys() if
                                    "weight_contact" in parameter_name]

            if (len(weight_names_bg) > 0):

                # first weight is set to 1 for all 400 ab
                name_weight_fixed = [parameter_name for parameter_name in weight_names_bg if "0" in parameter_name][0]
                # set weight1 to 1 because of overparametrizaton of softmax
                transformed_parameters[name_weight_fixed] = [1.0]

                # softmax
                for weight in [parameter_name for parameter_name in weight_names_bg if "0" not in parameter_name]:
                    transformed_parameters[weight] = np.log(parameters_in[weight][0] * np.exp(1)
                                                            / parameters_in[name_weight_fixed]).tolist()

            if (len(weight_names_contact) > 0):
                # first weight is set to 1 for all 400 ab
                name_weight_fixed = \
                [parameter_name for parameter_name in weight_names_contact if "0" in parameter_name][0]
                transformed_parameters[name_weight_fixed] = [
                    1.0]  # set weight1 to 1 because of overparametrizaton of softmax

                # softmax
                for weight in [parameter_name for parameter_name in weight_names_contact if "0" not in parameter_name]:
                    transformed_parameters[weight] = np.log(
                        parameters_in[weight][0] * np.exp(1) / parameters_in[name_weight_fixed]).tolist()

        # mean values will be transformed into log space
        if (mean):
            mean_columns = [parameter_name for parameter_name in parameters_in.keys() if "mu" in parameter_name]

            for mean in mean_columns:
                transformed_parameters[mean] = np.log(parameters_in[mean]).tolist()

        # precision will be transformed into log precision
        if (prec):
            prec_columns = [parameter_name for parameter_name in parameters_in.keys() if "prec" in parameter_name]

            for name in prec_columns:

                parameter = parameters_in[name]

                if (len(parameter) == 1):
                    transformed_parameters[name] = [np.log(parameter[0])]
                elif (len(parameter) == 400):
                    transformed_parameters[name] = np.log(parameter).tolist()  # row-wise
                else:
                    print('specified covariance matrix is neither diagonal nor isotrope. Do nothing. ')

    else:
        # weights will be backtransformed from softmax representation
        if (weights):

            weight_names_bg = [parameter_name for parameter_name in parameters_in.keys() if
                               "weight_bg" in parameter_name]
            weight_names_contact = [parameter_name for parameter_name in parameters_in.keys() if
                                    "weight_contact" in parameter_name]

            if (len(weight_names_bg) > 0):
                sum_weights = 0
                for weight in weight_names_bg:
                    transformed_parameters[weight] = np.exp(parameters_in[weight][0]).tolist()
                    sum_weights += transformed_parameters[weight]

                for weight in weight_names_bg:
                    transformed_parameters[weight] = [transformed_parameters[weight] / sum_weights]

            if (len(weight_names_contact) > 0):
                sum_weights = 0
                for weight in weight_names_contact:
                    transformed_parameters[weight] = np.exp(parameters_in[weight][0]).tolist()
                    sum_weights += transformed_parameters[weight]

                for weight in weight_names_contact:
                    transformed_parameters[weight] = [transformed_parameters[weight] / sum_weights]

        # log mean values will be transformed back via exp
        if (mean):
            mean_columns = [parameter_name for parameter_name in parameters_in.keys() if "mu" in parameter_name]

            for mean in mean_columns:
                transformed_parameters[mean] = np.exp(parameters_in[mean]).tolist()

        # log prec values will be transformed back via exp
        if (prec):
            prec_columns = [parameter_name for parameter_name in parameters_in.keys() if "prec" in parameter_name]

            for name in prec_columns:

                parameter = parameters_in[name]

                if (len(parameter) == 1):
                    transformed_parameters[name] = [np.exp(parameter[0])]
                elif (len(parameter) == 400):
                    transformed_parameters[name] = np.exp(parameter).tolist()
                else:
                    print('specified covariance matrix is neither diagonal nor full. Do nothing. ')

    return transformed_parameters

