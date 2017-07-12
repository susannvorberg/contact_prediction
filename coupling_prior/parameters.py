#!/usr/bin/env python



import numpy as np
import json


class Parameters():
    """
    Specify the Parameters for the Mixture Model of the Coupling Prior

    """

    def __init__(self, parameter_dir):

        #parameters
        self.parameters_structured          = {}
        self.parameters_linear              = np.array([])

        #path to output files
        self.parameter_dir          = str(parameter_dir)
        self.parameter_file         = self.parameter_dir + "/parameters"


        #parameter settings
        self.sigma                          = 'diagonal'
        self.nr_components                  = 4
        self.fixed_parameters               = ['weight_bg_0', 'weight_contact_0', 'mu_0', 'prec_0']
        self.prec_wrt_L                     = False

    def __repr__(self):


        str = "\nParameters of Gaussian Mixture for Coupling Prior:\n"

        str += "\nPath to ouput files: \n"
        for param in ["self.parameter_file"]:
            str += "{0:{1}} {2}\n".format(param.split(".")[1], '36', eval(param))


        str += "\nParameter settings: \n"
        for param in ["self.sigma",
                      "self.nr_components",
                      "self.fixed_parameters",
                      "self.prec_wrt_L"  ]:
            str += "{0:{1}} {2}\n".format(param.split(".")[1], '36', eval(param))


        return(str)


    def set_nr_components(self, nr_components):
        if nr_components > 0:
            self.nr_components = int(nr_components)

    def set_sigma(self, sigma, prec_wrt_L):

        if sigma in ['diagonal', 'isotrope', 'full']:
            self.sigma = sigma
        else:
            print("Sigma must be one of: ['diagonal', 'isotrope', 'full'].")
            self.sigma = 'diagonal'

        self.prec_wrt_L = prec_wrt_L

    def set_fixed_parameters(self, fixed_parameters):

        if not isinstance(fixed_parameters, list):
            print("fixed_parameters must be a list of parameter names.")
            self.fixed_parameters = ['weight_bg_0', 'weight_contact_0']
            return

        self.fixed_parameters = fixed_parameters

    def _initialise_parameters_default(self):

        #initialise parameters with a strong component at zero
        #with equal weights for all components
        initial_means       = [0]       + [0] * (self.nr_components-1)
        initial_weights     = [1.0 / self.nr_components] *  self.nr_components
        initial_precision   = [1/0.0005]  + [1/0.05] * (self.nr_components-1) #precision  = 1/variance
        if self.prec_wrt_L:
            initial_precision   = [1/1.0] + [1/0.8]* (self.nr_components-1) #precision  = 1/variance * L


        parameters = {}

        for component in range(len(initial_weights)):

            parameters['weight_contact_' + str(component)] = [initial_weights[component]]  ##* 400
            parameters['weight_bg_' + str(component)] = [initial_weights[component]]  ##* 400


            if 'mu_' + str(component) in self.fixed_parameters:
                parameters['mu_0'] = [initial_means[0]] * 400
            else:
                parameters['mu_' + str(component)] = np.random.normal(
                    loc=initial_means[component],
                    scale=0.05,
                    size=400
                ).tolist()

            if ('prec_' + str(component) in self.fixed_parameters) or self.sigma == 'isotrope':
                parameters['prec_' + str(component)] = [initial_precision[component]] * 400
            else:
                parameters['prec_' + str(component)] = np.random.normal(
                    loc=initial_precision[component],
                    scale=initial_precision[component] / 10,
                    size=400
                ).tolist()

        return parameters

    def initialise_parameters(self, nr_components=None, sigma=None, prec_wrt_L=None, fixed_parameters=None):

        if nr_components is not None:
            self.set_nr_components(nr_components)

        if sigma is not None and prec_wrt_L is not None:
            self.set_sigma(sigma, prec_wrt_L)

        if fixed_parameters is not None:
            self.set_fixed_parameters(fixed_parameters)


        print("\nInitialise parameters with default values using:")
        print("\tnr components: {0}".format(self.nr_components))
        print("\tprecision matrix is {0} and it's determined wrt to L: {1}".format(self.sigma, self.prec_wrt_L))
        print("\tfixed parameters that will not be optimized: {0}".format(self.fixed_parameters))

        self.parameters_structured =  self._initialise_parameters_default()

        #transform weights and precision
        self.parameters_structured = self.transform_parameters(weights=True, mean=False, prec=True, back=False)
        self.parameters_linear = self.structured_to_linear(self.parameters_structured)

        self.print_parameters()

    def print_parameters(self):

        print("\nParameters for Gaussian mixture model of coupling prior:")
        for key, value in self.parameters_structured.iteritems():
            if "weight" in key:
                print("{0}: \t {1}".format(key, value[0]))
        for key, value in self.parameters_structured.iteritems():
            if "mu" in key:
                print("{0}: \t\t\t min: {1} \t mean: {2} \t max: {3}".format(key,
                                                                             np.round(np.min(value), decimals=3),
                                                                             np.round(np.mean(value), decimals=3),
                                                                             np.round(np.max(value), decimals=3)))
        for key, value in self.parameters_structured.iteritems():
            if "prec" in key:
                print("{0}: \t\t min: {1} \t mean: {2} \t max: {3}".format(key,
                                                                           np.round(np.min(value), decimals=3),
                                                                           np.round(np.mean(value), decimals=3),
                                                                           np.round(np.max(value), decimals=3)))

    def write_parameters(self, parameters):

        with open(self.parameter_file, 'w') as fp:
            json.dump(parameters, fp)

    def linear_to_structured(self, parameters_linear):

        sorted_keys = [key for key in sorted(self.parameters_structured.keys()) if key not in self.fixed_parameters]

        parameters_structured = {}
        for position, key in enumerate(sorted_keys):
            if 'prec' in key and self.sigma == 'isotrope':
                parameters_structured[key]  = [parameters_linear[0]] * 400
                parameters_linear           = parameters_linear[1:]
            else:
                parameters_structured[key]  = parameters_linear[:len(self.parameters_structured[key])].tolist()
                parameters_linear           = parameters_linear[len(self.parameters_structured[key]):]

        for key in self.fixed_parameters:
            parameters_structured[key] = self.parameters_structured[key]

        return parameters_structured

    def structured_to_linear(self, parameters_structured):

        parameters_linear = []
        for key, values in sorted(parameters_structured.iteritems()):
            if key not in self.fixed_parameters:
                if 'prec' in key and self.sigma == 'isotrope':
                    parameters_linear.extend([values[0]])
                else:
                    parameters_linear.extend(values)

        return np.array(parameters_linear)

    def get_parameters_linear(self):
        return(self.parameters_linear)

    def get_parameters_structured(self):
        return(self.parameters_structured)

    def get_settings(self):
        settings = {}
        settings['nr_components'] = self.nr_components
        settings['sigma'] = self.sigma
        settings['prec_wrt_L'] = self.prec_wrt_L
        settings['fixed_parameters'] = self.fixed_parameters
        settings['parameter_file'] = self.parameter_file

        return settings

    def transform_parameters(self, weights=True, mean=False, prec=True, back=False):
        """
        parameters: Pandas dataframe with column names
        weights:    will be transformed with softmax
        mean:       will be transformed with exp
        prec:       will be transformed with exp
        back:       will regain original parameters
        """

        transformed_parameters = self.parameters_structured.copy()

        if (not back):

            # weights will be represented as softmax
            if (weights):
                weight_names_bg = [parameter_name for parameter_name in self.parameters_structured.keys() if
                                   "weight_bg" in parameter_name]
                weight_names_contact = [parameter_name for parameter_name in self.parameters_structured.keys() if
                                        "weight_contact" in parameter_name]

                if (len(weight_names_bg) > 0):

                    # first weight is set to 1 for all 400 ab
                    name_weight_fixed = [parameter_name for parameter_name in weight_names_bg if "0" in parameter_name][
                        0]
                    # set weight1 to 1 because of overparametrizaton of softmax
                    transformed_parameters[name_weight_fixed] = [1.0]

                    # softmax
                    for weight in [parameter_name for parameter_name in weight_names_bg if "0" not in parameter_name]:
                        transformed_parameters[weight] = np.log(self.parameters_structured[weight][0] * np.exp(1)
                                                                / self.parameters_structured[name_weight_fixed]).tolist()

                if (len(weight_names_contact) > 0):
                    # first weight is set to 1 for all 400 ab
                    name_weight_fixed = \
                        [parameter_name for parameter_name in weight_names_contact if "0" in parameter_name][0]
                    transformed_parameters[name_weight_fixed] = [
                        1.0]  # set weight1 to 1 because of overparametrizaton of softmax

                    # softmax
                    for weight in [parameter_name for parameter_name in weight_names_contact if
                                   "0" not in parameter_name]:
                        transformed_parameters[weight] = np.log(
                            self.parameters_structured[weight][0] * np.exp(1) / self.parameters_structured[name_weight_fixed]).tolist()

            # mean values will be transformed into log space
            if (mean):
                mean_columns = [parameter_name for parameter_name in self.parameters_structured.keys() if "mu" in parameter_name]

                for mean in mean_columns:
                    transformed_parameters[mean] = np.log(self.parameters_structured[mean]).tolist()

            # precision will be transformed into log precision
            if (prec):
                prec_columns = [parameter_name for parameter_name in self.parameters_structured.keys() if "prec" in parameter_name]

                for name in prec_columns:
                    parameter = self.parameters_structured[name]
                    transformed_parameters[name] = np.log(parameter).tolist()  # row-wise


        else:
            # weights will be backtransformed from softmax representation
            if (weights):

                weight_names_bg = [parameter_name for parameter_name in self.parameters_structured.keys() if
                                   "weight_bg" in parameter_name]
                weight_names_contact = [parameter_name for parameter_name in self.parameters_structured.keys() if
                                        "weight_contact" in parameter_name]

                if (len(weight_names_bg) > 0):
                    sum_weights = 0
                    for weight in weight_names_bg:
                        transformed_parameters[weight] = np.exp(self.parameters_structured[weight][0]).tolist()
                        sum_weights += transformed_parameters[weight]

                    for weight in weight_names_bg:
                        transformed_parameters[weight] = [transformed_parameters[weight] / sum_weights]

                if (len(weight_names_contact) > 0):
                    sum_weights = 0
                    for weight in weight_names_contact:
                        transformed_parameters[weight] = np.exp(self.parameters_structured[weight][0]).tolist()
                        sum_weights += transformed_parameters[weight]

                    for weight in weight_names_contact:
                        transformed_parameters[weight] = [transformed_parameters[weight] / sum_weights]

            # log mean values will be transformed back via exp
            if (mean):
                mean_columns = [parameter_name for parameter_name in self.parameters_structured.keys() if "mu" in parameter_name]

                for mean in mean_columns:
                    transformed_parameters[mean] = np.exp(self.parameters_structured[mean]).tolist()

            # log prec values will be transformed back via exp
            if (prec):
                prec_columns = [parameter_name for parameter_name in self.parameters_structured.keys() if "prec" in parameter_name]

                for name in prec_columns:
                    parameter = self.parameters_structured[name]
                    transformed_parameters[name] = np.exp(parameter).tolist()

        return transformed_parameters
