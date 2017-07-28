#!/usr/bin/env python


import json
import time
import sys
import coupling_prior.ext.libll as libll
import numpy as np
from statsmodels.tools import numdiff
import raw
import plotting.plot_optimization_logfile as opt_plot
import utils_coupling_prior as coupling_prior_utils
import utils.io_utils as io
from coupling_prior.likelihood_protein import LikelihoodProtein
from coupling_prior.parameters import Parameters
from coupling_prior.regularization import Regularization
import copy

class LikelihoodFct():
    """
    Specify the Likelihood and Gradients

    """

    def __init__(self, plot_dir ):

        self.dataset = None
        self.training_data  = {}
        self.test_data      = {}
        self.evaluation_set = {}
        self.batch_nr = 1

        #parameters
        self.parameters          = None

        #regularization settings
        self.regularization = None

        #path to output files
        self.settings_file          = None
        self.optimization_log_file  = None
        self.plot_dir               = str(plot_dir)
        self.plot_name              = self.plot_dir + "/optimization_monitor_plot.html"


        # settings
        self.debug_mode = 0
        self.threads_per_protein = 1  # no parallelized by default
        self.hessian_pseudocount = 0  # only for debugging


    def __repr__(self):


        str = "\nLikelihood Function and Gradients:\n"

        str += "\nPath to ouput files: \n"
        for param in ["self.plot_name",
                      "self.settings_file",
                      "self.optimization_log_file"]:
            str += "{0:{1}} {2}\n".format(param.split(".")[1], '36', eval(param))


        str += "\nComputation settings: \n"
        str += "{0:{1}} {2}\n".format("debug_mode", '36', self.debug_mode)
        str += "{0:{1}} {2}\n".format("threads_per_protein", '36', self.threads_per_protein)
        str += "{0:{1}} {2}\n".format("hessian_pseudocount", '36', self.hessian_pseudocount)

        if self.dataset:
            str += self.dataset.__repr__()

        if self.parameters:
            str += self.parameters.__repr__()

        if self.regularization:
            str += self.regularization.__repr__()


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

        if 'prec_wrt_L' in settings:
            self.set_prec_wrt_L = settings['prec_wrt_L']

    def get_settings(self):

        settings = {}
        settings['debug_mode'] = self.debug_mode
        settings['threads_per_protein'] = self.threads_per_protein
        settings['hessian_pseudocount'] = self.hessian_pseudocount
        settings['plot_name'] = self.plot_name
        settings['optimization_log_file'] = self.optimization_log_file
        settings['settings_file'] = self.settings_file

        if self.dataset:
            settings.update(self.dataset.get_settings())

        if self.parameters:
            settings.update(self.parameters.get_settings())

        if self.regularization:
            settings.update(self.regularization.get_settings())


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

    def set_parameters(self, parameters):
        self.parameters = parameters

        self.settings_file = self.parameters.parameter_file + ".settings"
        self.optimization_log_file = self.parameters.parameter_file + ".log"

    def set_debug_mode(self, debug_mode):
        self.debug_mode = int(debug_mode)

    def set_nr_threads_per_protein(self, nr_threads):
        self.threads_per_protein = int(nr_threads)

    def set_regularizer(self, reg_coeff_mu, reg_coeff_diag_prec):
        self.regularization = Regularization(self.parameters, reg_coeff_mu, reg_coeff_diag_prec)

    def write_settings(self):

        if self.settings_file is not None:
            with open(self.settings_file, 'w') as fp:
                json.dump(self.get_settings(), fp)

    def compute_regularization(self, parameters_structured):

        if self.regularization is None:
            print("You need to set the regularizer first with  'set_regularizer'. ")
            return

        #update regularizer with current parameters
        self.regularization.set_parameters(parameters_structured)

        #compute regularization
        reg = self.regularization.get_regularization()

        #compute gradients of regularizer
        reg_gradients_struct = self.regularization.get_regularization_gradients()

        return reg, reg_gradients_struct

    def compute_f(self, parameter, parameter_name, parameter_index, LikelihoodProtein):


        parameters_structured = copy.deepcopy(self.parameters.parameters_structured)

        if ('prec' in parameter_name) and (self.parameters.sigma == 'isotrope'):
            parameters_structured[parameter_name] = [parameter[0]] * 400
        else:
            parameters_structured[parameter_name][parameter_index] = parameter[0]

        #compute function value at current parameter value
        LikelihoodProtein.set_parameters(parameters_structured)
        LikelihoodProtein.compute_f()
        f = LikelihoodProtein.get_f()

        # regularization
        reg, reg_gradients = self.compute_regularization(parameters_structured)
        f = reg

        return f

    def compute_f_py(self, parameter, braw, Nij, qij, protein, parameter_name, parameter_index):


        #compute function value at current parameter value
        LL = LikelihoodProtein(braw, Nij, qij)
        LL.set_pairs(
            self.training_data[protein]['residue_i'],
            self.training_data[protein]['residue_j'],
            self.training_data[protein]['contact']
        )

        # define parameters
        parameters = Parameters(self.parameters.parameter_dir)
        parameters.set_sigma(self.parameters.sigma, self.parameters.prec_wrt_L)
        parameters.set_nr_components(self.parameters.nr_components)
        parameters.set_fixed_parameters(self.parameters.fixed_parameters)
        parameters.initialise_parameters(seed=123, verbose=False)

        parameters_structured = copy.deepcopy(self.parameters.parameters_structured)
        if ('prec' in parameter_name) and (self.parameters.sigma == 'isotrope'):
            parameters_structured[parameter_name] = [parameter[0]] * 400
        else:
            parameters_structured[parameter_name][parameter_index] = parameter[0]
        parameters.parameters_structured=parameters_structured

        #will also be transformed!
        LL.set_parameters(parameters)

        LL.compute_f_df(compute_gradients=False)
        f = LL.get_f()

        return f

    def numerical_gradient(self, check_weights=True, check_mu=True, check_prec=True):

        #select data for checking the gradient
        protein = self.training_data.keys()[0]
        data_for_gradient_check = {protein: self.training_data[protein]}

        print("\nUse protein {0} for gradient checking".format(protein))
        print("with {0} contacts and {1} non-contacts\n".format(
            np.sum(data_for_gradient_check[protein]['contact']),
            len(data_for_gradient_check[protein]['contact']) - np.sum(data_for_gradient_check[protein]['contact']))
        )


        ## compute analytical gradients
        print("compute analytical gradient...\n")
        LL = libll.Likelihood_Dataset(
            data_for_gradient_check,
            self.parameters.parameters_structured,
            self.parameters.prec_wrt_L
        )
        LL.set_debug_mode(0)
        LL.set_threads_per_protein(self.threads_per_protein)
        LL.compute_f_df(self.hessian_pseudocount)
        gradients = LL.get_gradient_dict()
        if self.parameters.sigma == 'isotrope':
            for key in gradients.keys():
                if 'prec' in key:
                    gradients[key] = [np.sum(gradients[key])] * 400


        # regularization
        self.regularization.fixed_parameters = [] #so that we also get gradients for fixed parameters
        reg, reg_gradients = self.compute_regularization(self.parameters.parameters_structured)
        for key in reg_gradients.keys():
            print key
            gradients[key] = reg_gradients[key]


        # compute analytical gradients
        # print("compute analytical gradient...\n")
        # braw = raw.parse_msgpack(self.training_data[protein]['braw_file_path'])
        # Nij, qij = io.read_qij(self.training_data[protein]['qijabfilename'], braw.ncol)
        # LL = LikelihoodProtein(braw, Nij, qij)
        # LL.set_pairs(
        #     self.parameters.parameters_structured[protein]['residue_i'],
        #     self.parameters.parameters_structured[protein]['residue_j'],
        #     self.parameters.parameters_structured[protein]['contact']
        # )
        # LL.set_parameters(parameters)
        # LL.compute_f_df()
        # gradients = LL.get_gradients()

        if(check_weights):
            print("\ncheck gradients of weights parameters...\n")

            print("{0:>20} {1:>20} {2:>20} {3:>20}".format(
                "parameter", "numdiff", "analytical", "|difference|"))

            for comp in range(self.parameters.nr_components):

                parameter = 'weight_bg_'+str(comp)
                init_param=copy.deepcopy(self.parameters.parameters_structured[parameter][0])
                num_grad_weight = numdiff.approx_fprime(
                    [init_param],
                    self.compute_f,
                    args=(parameter, 0, LL)
                )
                diff = abs(gradients[parameter][0] - num_grad_weight)
                print("{0:>20} {1:>20} {2:>20} {3:>20}".format(
                    parameter, num_grad_weight, gradients[parameter][0], diff))


                parameter = 'weight_contact_' + str(comp)
                init_param=copy.deepcopy(self.parameters.parameters_structured[parameter][0])
                num_grad_weight = numdiff.approx_fprime(
                    [init_param],
                    self.compute_f,
                    args=(parameter, 0, LL)
                )
                diff = abs(gradients[parameter][0] - num_grad_weight)
                print("{0:>20} {1:>20} {2:>20} {3:>20}".format(
                    parameter, num_grad_weight, gradients[parameter][0], diff))


        if(check_mu):
            print("\ncheck gradients of mean parameters...\n")

            print("{0:>20} {1:>20} {2:>20} {3:>20} {4:>20}".format(
                "parameter", "index", "numdiff", "analytical", "|difference|"))

            parameter_indices = np.random.randint(400, size=3)

            for comp in range(self.parameters.nr_components):
                for parameter_index in parameter_indices.tolist():

                    parameter = 'mu_'+str(comp)
                    init_param = copy.deepcopy(self.parameters.parameters_structured[parameter][parameter_index])
                    num_grad_weight = numdiff.approx_fprime(
                        [init_param],
                        self.compute_f,
                        args=(parameter, parameter_index, LL)
                    )
                    diff = abs(gradients[parameter][parameter_index] - num_grad_weight)
                    print("{0:>20} {1:>20} {2:>20} {3:>20} {4:>20}".format(
                        parameter, parameter_index, num_grad_weight, gradients[parameter][parameter_index], diff))


        if(check_prec):
            print("\ncheck gradients of precision parameters...\n")

            print("{0:>20} {1:>20} {2:>20} {3:>20} {4:>20}".format(
                "parameter", "index", "numdiff", "analytical", "|difference|"))

            parameter_indices = np.random.randint(400, size=3)

            for comp in range(self.parameters.nr_components):
                for parameter_index in parameter_indices.tolist():

                    parameter = 'prec_'+str(comp)
                    init_param = copy.deepcopy(self.parameters.parameters_structured[parameter][parameter_index])
                    # num_grad_weight = numdiff.approx_fprime(
                    #     [init_param],
                    #     self.compute_f_py,
                    #     args=(braw, Nij, qij, protein, parameter, parameter_index)
                    # )
                    num_grad_weight = numdiff.approx_fprime(
                        [init_param],
                        self.compute_f,
                        args=(parameter, parameter_index, LL)
                    )
                    diff = abs(gradients[parameter][parameter_index] - num_grad_weight)
                    print("{0:>20} {1:>20} {2:>20} {3:>20} {4:>20}".format(
                        parameter, parameter_index, num_grad_weight, gradients[parameter][parameter_index], diff))

    @staticmethod
    def compute_neg_log_likelihood_protein(braw_file, qij_file, residues_i, residues_j, parameters, contact=1):

        assert(len(residues_i) == len(residues_j)), "residues_i and residues_j are not of same length!"

        braw = raw.parse_msgpack(braw_file)
        Nij, qij = io.read_qij(qij_file, braw.ncol)

        lik_protein = LikelihoodProtein(braw, Nij, qij)

        lik_protein.set_pairs(
            residues_i,
            residues_j,
            [contact] * len(residues_i)
        )

        lik_protein.set_parameters(parameters)

        lik_protein.compute_f_df(compute_gradients=False)

        f_protein = lik_protein.get_f_pairwise()

        return f_protein

    def f_df_protein(self, protein, parameters_structured):

        print(protein)

        braw = raw.parse_msgpack(self.training_data[protein]['braw_file_path'])
        Nij, qij = io.read_qij(self.training_data[protein]['qijabfilename'], braw.ncol)

        lik_protein = LikelihoodProtein(braw, Nij, qij)

        lik_protein.set_pairs(
            self.training_data[protein]['residue_i'],
            self.training_data[protein]['residue_j'],
            self.training_data[protein]['contact']
        )

        #create new Parameter instance
        parameters = Parameters(self.parameters.parameter_dir)
        parameters.set_sigma(self.parameters.sigma, self.parameters.prec_wrt_L)
        parameters.set_nr_components(self.parameters.nr_components)
        parameters.set_fixed_parameters(self.parameters.fixed_parameters)
        parameters.set_parameters_structured(parameters_structured)
        lik_protein.set_parameters(parameters)

        lik_protein.compute_f_df()

        f_protein = lik_protein.get_f()
        grad_protein = lik_protein.get_gradients()

        return f_protein, grad_protein

    def f_df_py(self, parameters_linear):

        f = 0
        g = {}

        parameters_linear = np.array(parameters_linear)
        parameters_structured = self.parameters.linear_to_structured(parameters_linear)

        #initialise gradient dict with zeros
        for parameter in self.parameters.parameters_structured.keys():
            if parameter not in self.parameters.fixed_parameters:
                g[parameter] = [0] * len(self.parameters.parameters_structured[parameter])

        t = time.time()
        for protein in self.training_data.keys():
            f_protein, grad_protein = self.f_df_protein(protein, parameters_structured.copy())
            f += f_protein
            for parameter, val in grad_protein.iteritems():
                g[parameter] += val
        timestamp = time.time() - t

        self.print_status(f, g, timestamp)

        ##### update log file
        log_df = self.update_optimization_logfile(f, g, timestamp)

        parameters_transformed_back = self.parameters.transform_parameters(weights=True, mean=False, prec=True, back=True)

        ##### save parameters and settings
        self.parameters.write_parameters(parameters_transformed_back)
        self.write_settings()

        ##### Plot the evaluation plot
        # opt_plot.plot_evaluation(parameters_transformed_back,
        #                          log_df,
        #                          self.get_settings(),
        #                          self.evaluation_set,
        #                          self.plot_name
        #                          )

        return(f, self.parameters.structured_to_linear(g))

    def f_df(self, parameters_linear):

        self.parameters.parameters_linear = np.array(parameters_linear)
        self.parameters.parameters_structured = self.parameters.linear_to_structured(parameters_linear)

        LL = libll.Likelihood_Dataset(
            self.training_data,
            self.parameters.parameters_structured,
            self.parameters.prec_wrt_L
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
        reg, reg_gradients = self.compute_regularization(self.parameters.parameters_structured)
        f -= reg
        for key in reg_gradients.keys():
            gradients[key] =  list(np.array(gradients[key]) - np.array(reg_gradients[key]))

        #handle isotrope case: gradients of all 400 dimenstions need to be summed
        if self.parameters.sigma == 'isotrope':
            for key in gradients.keys():
                if 'prec' in key:
                    gradients[key] = [np.sum(gradients[key])] * 400

        self.print_status(f+reg, reg,  gradients, timestamp)

        ##### update log file
        log_df = self.update_optimization_logfile(f, gradients, timestamp)


        parameters_transformed_back = self.parameters.transform_parameters(weights=True, mean=False, prec=True, back=True)

        ##### save parameters and settings
        self.parameters.write_parameters(parameters_transformed_back)
        self.write_settings()

        ##### Plot the evaluation plot
        opt_plot.plot_evaluation(parameters_transformed_back,
                                 log_df,
                                 self.get_settings(),
                                 self.evaluation_set,
                                 self.plot_name
                                 )

        return(f, self.parameters.structured_to_linear(gradients))

    @staticmethod
    def print_status(f, reg, gradients_dict, timestamp):
        gradient_norm_weight = np.linalg.norm(
            [x for key, value in gradients_dict.iteritems() if 'weight' in key for x in value])
        gradient_norm_mu = np.linalg.norm([x for key, value in gradients_dict.iteritems() if 'mu' in key for x in value])
        gradient_norm_prec = np.linalg.norm(
            [x for key, value in gradients_dict.iteritems() if 'prec' in key for x in value])

        header_tokens = [("f", 24),
                         ("reg", 24),
                         ("gradient_norm_weight", 24),
                         ("gradient_norm_mu", 24),
                         ("gradient_norm_prec", 24),
                         ("timestamp", 24)]
        headerline = (" ".join("{0:>{1}s}".format(ht, hw) for ht, hw in header_tokens))
        print("\n" + headerline)


        print("{0:>24} {1:>24} {2:>24} {3:>24} {4:>24} {5:>24}\n".format(
            f, reg, gradient_norm_weight, gradient_norm_mu, gradient_norm_prec, timestamp
        ))

        sys.stdout.flush()

    def update_optimization_logfile(self, f, gradients, timestamp):

        log_df = coupling_prior_utils.read_optimization_log_file(self.optimization_log_file)

        f_crossval = np.nan
        if len(log_df) % 10 == 0:
            optim_cpp = libll.Likelihood_Dataset(
                self.test_data,
                self.parameters.parameters_structured,
                self.parameters.prec_wrt_L
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

