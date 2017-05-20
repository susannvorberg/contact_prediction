#!/usr/bin/env python


import os
import numpy as np
import time
import json
import pandas as pd
from statsmodels.tools import numdiff

import utils_coupling_prior as coupling_prior_utils
import plotting.plot_optimization_logfile as opt_plot
import coupling_prior.ext.libreg as libreg
import coupling_prior.ext.libll as libll

class LikelihoodFct():
    """
    Specify the Likelihood and Gradients

    """

    def __init__(self, parameter_dir, plot_dir ):

        self.dataset = None
        self.training_data  = {}
        self.test_data      = {}
        self.evaluation_set = {}
        self.batch_nr = 1

        #parameters
        self.parameters_structured          = {}
        self.parameters_linear              = np.array([])

        #path to output files
        self.parameter_dir          = str(parameter_dir)
        self.parameter_file         = self.parameter_dir + "/parameters"
        self.settings_file          = self.parameter_file + ".settings"
        self.optimization_log_file  = self.parameter_file + ".log"
        self.plot_dir               = str(plot_dir)
        self.plot_name              = self.plot_dir + "/optimization_monitor_plot.html"

        #parameter settings
        self.sigma                          = 'diagonal'
        self.nr_components                  = 4
        self.regularizer_mu                 = 0 #set to 0 to ignore regularization
        self.regularizer_diagonal_precMat   = 0 #set to 0 to ignore regularization
        self.hessian_pseudocount            = 0 #only for debugging
        self.fixed_parameters               = ['weight_bg_0', 'weight_contact_0', 'mu_0', 'prec_0']

        # settings
        self.debug_mode = 0
        self.threads_per_protein = 1  # no parallelized by default


    def __repr__(self):


        str = "\nLikelihood Function and Gradients:\n"

        str += "\nPath to ouput files: \n"
        for param in ["self.parameter_file",
                      "self.plot_name",
                      "self.settings_file",
                      "self.optimization_log_file"]:
            str += "{0:{1}} {2}\n".format(param.split(".")[1], '36', eval(param))


        str += "\nParameter settings: \n"
        for param in ["self.sigma",
                      "self.nr_components",
                      "self.regularizer_mu",
                      "self.regularizer_diagonal_precMat",
                      "self.hessian_pseudocount",
                      "self.fixed_parameters"]:
            str += "{0:{1}} {2}\n".format(param.split(".")[1], '36', eval(param))


        str += "\nComputation settings: \n"
        str += "{0:{1}} {2}\n".format("debug_mode", '36', self.debug_mode)
        str += "{0:{1}} {2}\n".format("threads_per_protein", '36', self.threads_per_protein)



        return(str)

    def read_settings(self, settings_file ):

        settings = coupling_prior_utils.read_settings_file(settings_file)

        if 'nr_components' in settings:
            self.nr_components = settings['nr_components']

        if 'regularizer_mu' in settings:
            self.regularizer_mu = settings['regularizer_mu']

        if 'regularizer_diagonal_covMat' in settings:
            self.regularizer_diagonal_covMat = settings['regularizer_diagonal_covMat']

        if 'fixed_parameters' in settings:
            self.fixed_parameters = settings['fixed_parameters']

        if 'sigma' in settings:
            self.sigma = settings['sigma']

        if 'threads_per_protein' in settings:
            self.threads_per_protein = settings['threads_per_protein']

    def get_settings(self):

        settings = {}
        settings['debug_mode'] = self.debug_mode
        settings['threads_per_protein'] = self.threads_per_protein
        settings['fixed_parameters'] = self.fixed_parameters
        settings['regularizer_diagonal_precMat'] = self.regularizer_diagonal_precMat
        settings['regularizer_mu'] = self.regularizer_mu
        settings['nr_components'] = self.nr_components
        settings['sigma'] = self.sigma

        settings['plot_name'] = self.plot_name
        settings['optimization_log_file'] = self.optimization_log_file
        settings['settings_file'] = self.settings_file
        settings['parameter_file'] = self.parameter_file

        if self.dataset is not None:
            settings.update(self.dataset.get_settings())

        return(settings)

    def set_data(self, dataset, batch_nr=1):

        self.dataset = dataset
        self.training_data = self.dataset.get_training_data()
        self.test_data = self.dataset.get_test_data()

        couplings_contacts, couplings_noncontacts = self.dataset.get_decoy_set()
        self.evaluation_set['contact'] = np.array(couplings_contacts).transpose()
        self.evaluation_set['bg'] = np.array(couplings_noncontacts).transpose()

        #if we are dealing with stochastic optimization over mini batches
        self.batch_nr=batch_nr

    def set_debug_mode(self, debug_mode):
        self.debug_mode = int(debug_mode)

    def set_nr_threads_per_protein(self, nr_threads):
        self.threads_per_protein = int(nr_threads)

    def set_nr_components(self, nr_components):
        self.nr_components = int(nr_components)

    def set_regularizer_mu(self, reg_coeff_mu):
        self.regularizer_mu = float(reg_coeff_mu)

    def set_regularizer_diagonal_precMat(self, reg_coeff_diagPrec):
        self.regularizer_diagonal_precMat = float(reg_coeff_diagPrec)

    def set_sigma(self, sigma):
        self.sigma = sigma

        if sigma not in ['diagonal', 'isotrope', 'full']:
            print("Sigma must be one of: ['diagonal', 'isotrope', 'full'].")
            self.sigma = 'diagonal'

    def set_fixed_parameters(self, fixed_parameters):

        if not isinstance(fixed_parameters, list):
            print("fixed_parameters must be a list of parameter names.")
            self.fixed_parameters = ['weight_bg_0', 'weight_contact_0']
            return

        self.fixed_parameters = fixed_parameters


    def initialise_parameters(self, parameter_file=None,  nr_components=None):

        if parameter_file is not None:
            self.parameter_file = parameter_file
        if nr_components is not None:
            self.nr_components = nr_components


        if parameter_file is not None and os.path.exists(parameter_file):

            parameters  = coupling_prior_utils.read_parameter_file(parameter_file)
            self.read_settings(parameter_file + ".settings")
        else:

            #initialise parameters with a strong component at zero
            initial_means       = [0]       + [0] * (self.nr_components-1)
            initial_precision   = [1/0.0005]  + [1/0.05] * (self.nr_components-1) #precision  = 1/variance
            initial_weights     = [1.0 / self.nr_components] *  self.nr_components

            parameters = coupling_prior_utils.init_by_default(initial_weights, initial_means, initial_precision, self.sigma, self.fixed_parameters)

        self.settings_file = self.parameter_file + ".settings"
        self.optimization_log_file = self.parameter_file + ".log"

        parameters = coupling_prior_utils.transform_parameters(parameters, weights=True, mean=False, prec=True, back=False)

        self.parameters_structured = parameters
        self.parameters_linear = self.structured_to_linear(parameters)

        self.print_parameters()

    def print_parameters(self):

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

    def write_parameters_and_settings(self, parameters):

        with open(self.parameter_file, 'w') as fp:
            json.dump(parameters, fp)

        with open(self.settings_file, 'w') as fp:
            json.dump(self.get_settings(), fp)


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

    def regularizer(self):

        reg = 0
        reg_gradients_struct = {}

        if self.regularizer_mu != 0 or self.regularizer_diagonal_precMat != 0:

            regularizer = libreg.Regularizer(self.parameters_structured,
                                          self.regularizer_mu,
                                          self.regularizer_diagonal_precMat,
                                          self.debug_mode
                                          )
            regularizer.print_regularization_parameters()

            # Add regularization to likelihood (not for fixed component 0!!!)
            if self.regularizer_mu != 0:
                reg += regularizer.regularizer_mu()
            if self.regularizer_diagonal_precMat != 0:
                reg += regularizer.regularizer_diagPrecMat()

            # Add gradient of regularizer
            for parameter in sorted(self.parameters_structured):
                comp = int(parameter.split("_")[-1])  # 0 or 1
                if('mu' in parameter) and self.regularizer_mu != 0:
                    reg_gradients_struct[parameter] = np.array(regularizer.gradient_mu_comp_reg(comp))
                if('prec' in parameter) and self.regularizer_diagonal_precMat != 0:
                    reg_gradients_struct[parameter] =  np.diag(regularizer.gradient_diagprecMat_comp_reg(comp))


        return(reg, reg_gradients_struct)


    def compute_f(self, parameter, parameter_name, parameters_structured, parameter_index, LL):

        #parameter value to compute finite differences
        parameters_structured[parameter_name][parameter_index] = parameter[0]

        #compute function value at current parameter value
        LL.set_parameters(parameters_structured)
        LL.compute_f()
        f = LL.get_f()

        return f


    def compute_grad_approx(self, parameter_name, parameters_structured, parameter_index, LL):

        start_parameter_value = parameters_structured[parameter_name][parameter_index]

        numerical_gradient = float(
            numdiff.approx_fprime(
                [start_parameter_value],
                self.compute_f,
                args=(parameter_name, parameters_structured.copy(), parameter_index, LL)
            )
        )

        return numerical_gradient

    def numerical_gradient(self, check_weights=True, check_mu=True, check_prec=True):


        print("\nDataset used for gradient checking \nnumber of contacts: {0}, number of non-contacts: {1}\n".format(self.dataset.nr_pairs_contact, self.dataset.nr_pairs_noncontact))


        #compute analytical gradients
        print("compute analytical gradient...\n")
        LL = libll.Likelihood_Dataset(
            self.training_data,
            self.parameters_structured
        )
        LL.set_debug_mode(0)
        LL.set_threads_per_protein(self.threads_per_protein)
        LL.compute_f_df(self.hessian_pseudocount)
        gradients = LL.get_gradient_dict()




        if(check_weights):
            print("\ncheck gradients of weights parameters...\n")

            print("{0:>20} {1:>20} {2:>20} {3:>20}".format(
                "parameter", "numdiff", "analytical", "|difference|"))

            for comp in range(self.nr_components):
                parameter = 'weight_bg_'+str(comp)
                num_grad_weight= self.compute_grad_approx(parameter, self.parameters_structured, 0, LL)
                diff = abs(gradients[parameter][0] - num_grad_weight)
                print("{0:>20} {1:>20} {2:>20} {3:>20}".format(
                    parameter, num_grad_weight, gradients[parameter][0], diff))

                parameter = 'weight_contact_' + str(comp)
                num_grad_weight = self.compute_grad_approx(parameter, self.parameters_structured, 0, LL)
                diff = abs(gradients[parameter][0] - num_grad_weight)
                print("{0:>20} {1:>20} {2:>20} {3:>20}".format(
                    parameter, num_grad_weight, gradients[parameter][0], diff))


        if(check_mu):
            print("\ncheck gradients of mean parameters...\n")

            print("{0:>20} {1:>20} {2:>20} {3:>20} {4:>20}".format(
                "parameter", "index", "numdiff", "analytical", "|difference|"))

            parameter_indices = np.random.randint(400, size=20)

            for comp in range(self.nr_components):
                for parameter_index in parameter_indices:

                    parameter = 'mu_'+str(comp)
                    num_grad_weight= self.compute_grad_approx(parameter, self.parameters_structured, parameter_index, LL)
                    diff = abs(gradients[parameter][parameter_index] - num_grad_weight)
                    print("{0:>20} {1:>20} {2:>20} {3:>20} {4:>20}".format(
                        parameter, parameter_index, num_grad_weight, gradients[parameter][parameter_index], diff))


        if(check_prec):
            print("\ncheck gradients of precision parameters...\n")

            print("{0:>20} {1:>20} {2:>20} {3:>20} {4:>20}".format(
                "parameter", "index", "numdiff", "analytical", "|difference|"))

            parameter_indices = np.random.randint(400, size=20)

            for comp in range(self.nr_components):
                for parameter_index in parameter_indices:

                    parameter = 'prec_'+str(comp)
                    num_grad_weight= self.compute_grad_approx(parameter, self.parameters_structured, parameter_index, LL)
                    diff = abs(gradients[parameter][parameter_index] - num_grad_weight)
                    print("{0:>20} {1:>20} {2:>20} {3:>20} {4:>20}".format(
                        parameter, parameter_index, num_grad_weight, gradients[parameter][parameter_index], diff))


    def f_df(self, parameters_linear):

        self.parameters_linear = np.array(parameters_linear)
        self.parameters_structured = self.linear_to_structured(parameters_linear)

        LL = libll.Likelihood_Dataset(
            self.training_data,
            self.parameters_structured
        )

        LL.set_debug_mode(self.debug_mode)
        LL.set_threads_per_protein(self.threads_per_protein)

        # compute likelihood and gradients
        t = time.time()
        LL.compute_f_df(self.hessian_pseudocount)
        timestamp = time.time() - t

        # retrieve f and df
        f = LL.get_f()
        gradients = LL.get_gradient_dict()


        # regularization
        reg, reg_gradients = self.regularizer()
        f -= reg
        for key in reg_gradients.keys():
            gradients[key] -= reg_gradients[key]

        #handle isotrope case: gradients of all 400 dimenstions need to be summed
        if self.sigma == 'isotrope':
            for key in gradients.keys():
                if 'prec' in key:
                    gradients[key] = [np.sum(gradients[key])] * 400

        self.print_status(f, gradients, timestamp)

        ##### update log file
        log_df = self.update_optimization_logfile(f, gradients, timestamp)


        parameters_transformed_back = coupling_prior_utils.transform_parameters(
            self.parameters_structured, weights=True, mean=False, prec=True, back=True)

        ##### save parameters and settings
        self.write_parameters_and_settings(parameters_transformed_back)

        ##### Plot the evaluation plot
        opt_plot.plot_evaluation(parameters_transformed_back,
                                 log_df,
                                 self.get_settings(),
                                 self.evaluation_set,
                                 self.plot_name
                                 )

        return(f, self.structured_to_linear(gradients))


    @staticmethod
    def print_status(f, gradients_dict, timestamp):
        gradient_norm_weight = np.linalg.norm(
            [x for key, value in gradients_dict.iteritems() if 'weight' in key for x in value])
        gradient_norm_mu = np.linalg.norm([x for key, value in gradients_dict.iteritems() if 'mu' in key for x in value])
        gradient_norm_prec = np.linalg.norm(
            [x for key, value in gradients_dict.iteritems() if 'prec' in key for x in value])

        header_tokens = [("f", 24),
                         ("gradient_norm_weight", 24),
                         ("gradient_norm_mu", 24),
                         ("gradient_norm_prec", 24),
                         ("timestamp", 24)]
        headerline = (" ".join("{0:>{1}s}".format(ht, hw) for ht, hw in header_tokens))
        print("\n" + headerline)


        print("{0:>24} {1:>24} {2:>24} {3:>24} {4:>24} \n".format(
            f, gradient_norm_weight, gradient_norm_mu, gradient_norm_prec, timestamp
        ))

    def update_optimization_logfile(self, f, gradients, timestamp):

        log_df = coupling_prior_utils.read_optimization_log_file(self.optimization_log_file)

        f_crossval = np.nan
        if len(log_df) % 10 == 0:
            optim_cpp = libll.Likelihood_Dataset(
                self.test_data,
                self.parameters_structured
            )
            optim_cpp.set_debug_mode(1)
            optim_cpp.set_threads_per_protein(self.threads_per_protein)
            optim_cpp.compute_f()  # compute likelihood
            f_crossval = optim_cpp.get_f()  # get likelihood

        # compute GRADIENT NORMS OF ALL parameters
        gradient_norm_weight = np.linalg.norm([x for key, value in gradients.iteritems() if 'weight' in key for x in value])
        gradient_norm_mu = np.linalg.norm([x for key, value in gradients.iteritems() if 'mu' in key for x in value])
        gradient_norm_prec = np.linalg.norm([x for key, value in gradients.iteritems() if 'prec' in key for x in value])


        ##### Update log file
        new_log_entry = {
            'pass': 1,
            'step': len(log_df) + 1,
            'timestamp': timestamp,
            'negLL': f,
            'negLL_crossval': f_crossval,
            'gradient_norm_weights':gradient_norm_weight,
            'gradient_norm_means': gradient_norm_mu,
            'gradient_norm_prec':gradient_norm_prec
        }
        new_log_entry.update(gradients)
        log_df = log_df.append(new_log_entry, ignore_index=True)

        json_str = log_df.to_json(orient='split')
        with open(self.optimization_log_file, 'w') as outfile:
            json.dump(json_str, outfile)

        return log_df

