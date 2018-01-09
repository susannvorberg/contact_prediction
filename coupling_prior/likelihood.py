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
from  sklearn.neighbors import KernelDensity
from utils.io_utils import AB, AB_INDICES
import copy
import pathos.multiprocessing as mp
from functools import partial




class LikelihoodFct():
    """
    Specify the Likelihood and Gradients

    """

    def __init__(self, plot_dir ):

        self.dataset = None
        self.training_data  = {}
        self.test_data      = {}
        self.batch_nr = 1

        #plot options
        self.evaluation_set = {}
        self.size_evaluationset = 2000
        self.evaluation_set_kde ={}

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
        self.python_parallel = False
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

        #if we are dealing with stochastic optimization over mini batches
        self.batch_nr=batch_nr

    def initialize_optimization_plot(self):

        if self.dataset is None:
            print("Set data first before initializing plot options!")
            return


        if self.parameters is None:
            print("Parameter needs to be set!")
            return

        ab_list = [
            AB_INDICES['A-A'],
            AB_INDICES['C-C'],
            AB_INDICES['E-R'],
            AB_INDICES['R-E'],
            AB_INDICES['K-E'],
            AB_INDICES['E-E'],
            AB_INDICES['K-K'],
            AB_INDICES['K-R'],
            AB_INDICES['V-I'],
            AB_INDICES['I-L'],
            AB_INDICES['S-T'],
            AB_INDICES['S-S'],
            AB_INDICES['K-P'],
            AB_INDICES['N-N'],
            AB_INDICES['W-W'],
            AB_INDICES['G-F']
        ]

        couplings_contacts, couplings_noncontacts, avg_lambda_pair = self.dataset.get_decoy_set(size=self.size_evaluationset)
        self.evaluation_set['contact'] = np.array(couplings_contacts).transpose()
        self.evaluation_set['bg'] = np.array(couplings_noncontacts).transpose()

        bandwidth = 0.01
        self.evaluation_set_kde = {}
        self.evaluation_set_kde['x_grid'] = np.linspace(-0.5, 0.5, 500)
        self.evaluation_set_kde['contact'] = {}
        self.evaluation_set_kde['bg'] = {}

        # kernel density estimate for couplings wijab
        for ab in ab_list:
            kde_contact = KernelDensity(kernel='gaussian', bandwidth=bandwidth).fit(self.evaluation_set['contact'][ab].reshape(-1, 1))
            kde_bg = KernelDensity(kernel='gaussian', bandwidth=bandwidth).fit(self.evaluation_set['bg'][ab].reshape(-1, 1))

            ### add empirical distribution for example data points
            self.evaluation_set_kde['contact'][ab] = np.exp(kde_contact.score_samples(self.evaluation_set_kde['x_grid'].reshape(-1, 1)))
            self.evaluation_set_kde['bg'][ab] = np.exp(kde_bg.score_samples(self.evaluation_set_kde['x_grid'].reshape(-1, 1)))

        #sample points according to regularizer
        std_dev = np.sqrt(1.0/avg_lambda_pair)
        regularizer = np.random.normal(scale=std_dev, size=10000)
        kde_reg = KernelDensity(kernel='gaussian', bandwidth=0.1).fit(regularizer.reshape(-1, 1))
        self.evaluation_set_kde['regularizer'] = np.exp(kde_reg.score_samples(self.evaluation_set_kde['x_grid'].reshape(-1, 1)))

    def set_parameters(self, parameters):
        self.parameters = parameters

        self.settings_file = self.parameters.parameter_file + ".settings"
        self.optimization_log_file = self.parameters.parameter_file + ".log"

    def set_debug_mode(self, debug_mode, python_parallel=False):
        self.debug_mode = int(debug_mode)
        self.python_parallel = python_parallel

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

    def compute_f(self, parameter, parameter_name, parameter_index):

        #save original parameters
        parameters = copy.deepcopy(self.parameters)

        if ('prec' in parameter_name) and (parameters.sigma == 'isotrope'):
            self.parameters.parameters_structured[parameter_name] = [parameter[0]] * 400
        else:
            self.parameters.parameters_structured[parameter_name][parameter_index] = parameter[0]

        self.parameters_linear = self.parameters.structured_to_linear(self.parameters.parameters_structured)
        f, gradients = self.f_df(self.parameters_linear, plot=False, save=False)

        #restore original parameters
        self.parameters = parameters

        return f

    def compute_f_py(self, parameter, parameter_name, parameter_index):

        #save original parameters
        parameters = copy.deepcopy(self.parameters)

        if ('prec' in parameter_name) and (parameters.sigma == 'isotrope'):
            self.parameters.parameters_structured[parameter_name] = [parameter[0]] * 400
        else:
            self.parameters.parameters_structured[parameter_name][parameter_index] = parameter[0]

        self.parameters_linear = self.parameters.structured_to_linear(self.parameters.parameters_structured)
        f, _ = self.f_df_py(self.parameters_linear, plot=False, save=False, compute_gradients=False)

        #restore original parameters
        self.parameters = parameters

        return f

    def numerical_gradient(self, check_weights=True, check_mu=True, check_prec=True, use_py=False):

        #select data for checking the gradient
        protein = self.training_data.keys()[0]
        self.training_data = {protein: self.training_data[protein]}

        #so that we also get gradients for fixed parameters
        self.parameters.fixed_parameters = []
        self.regularization.fixed_parameters = []
        self.parameters.parameters_linear = self.parameters.structured_to_linear(self.parameters.parameters_structured)


        print("\nUse protein {0} for gradient checking".format(protein))
        print("with {0} contacts and {1} non-contacts\n".format(
            np.sum(self.training_data[protein]['contact']),
            len(self.training_data[protein]['contact']) - np.sum(self.training_data[protein]['contact']))
        )

        ##compute analytical gradients
        if use_py:

            ##compute analytical gradients
            print("compute analytical gradient with python function...\n")

            f, gradients_linear = self.f_df_py(self.parameters.parameters_linear, plot=False, save=False, compute_gradients=True)
            gradients = self.parameters.linear_to_structured(gradients_linear)
            f_value_function = self.compute_f_py

        else:
            ## compute analytical gradients
            print("compute analytical gradient with cpp function...\n")

            f, gradients_linear = self.f_df(self.parameters.parameters_linear, plot=False, save=False)
            gradients = self.parameters.linear_to_structured(gradients_linear)
            f_value_function = self.compute_f


        if(check_weights):
            print("\ncheck gradients of weights parameters...\n")

            print("{0:>20} {1:>20} {2:>20} {3:>20}".format(
                "parameter", "numdiff", "analytical", "|difference|"))

            for comp in range(self.parameters.nr_components):

                parameter = 'weight_bg_'+str(comp)
                if parameter in gradients:
                    init_param=copy.deepcopy(self.parameters.parameters_structured[parameter][0])
                    num_grad_weight = numdiff.approx_fprime(
                        [init_param],
                        f_value_function,
                        args=(parameter, 0)
                    )
                    diff = abs(gradients[parameter][0] - num_grad_weight)
                    print("{0:>20} {1:>20} {2:>20} {3:>20}".format(
                        parameter, num_grad_weight, gradients[parameter][0], diff))


                parameter = 'weight_contact_' + str(comp)
                if parameter in gradients:
                    init_param=copy.deepcopy(self.parameters.parameters_structured[parameter][0])
                    num_grad_weight = numdiff.approx_fprime(
                        [init_param],
                        f_value_function,
                        args=(parameter, 0)
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
                    if parameter in gradients:
                        init_param = copy.deepcopy(self.parameters.parameters_structured[parameter][parameter_index])
                        num_grad_weight = numdiff.approx_fprime(
                            [init_param],
                            f_value_function,
                            args=(parameter, parameter_index)
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
                    if parameter in gradients:
                        init_param = copy.deepcopy(self.parameters.parameters_structured[parameter][parameter_index])
                        # num_grad_weight = numdiff.approx_fprime(
                        #     [init_param],
                        #     self.compute_f_py,
                        #     args=(braw, Nij, qij, protein, parameter, parameter_index)
                        # )
                        num_grad_weight = numdiff.approx_fprime(
                            [init_param],
                            f_value_function,
                            args=(parameter, parameter_index)
                        )
                        diff = abs(gradients[parameter][parameter_index] - num_grad_weight)
                        print("{0:>20} {1:>20} {2:>20} {3:>20} {4:>20}".format(
                            parameter, parameter_index, num_grad_weight, gradients[parameter][parameter_index], diff))

    @staticmethod
    def f_df_py_protein_parallel(protein, parameters_transformed_back, prec_wrt_L, nr_components, fixed_parameters, status=True, compute_gradients=True):


        protein_name, protein_data = protein
        if status:
            print("Compute likelihood and gradients for {0}".format(protein_name))

        braw = raw.parse_msgpack(protein_data['braw_file_path'])
        Nij, qij = io.read_qij(protein_data['qijabfilename'], braw.ncol)
        L = braw.ncol

        #prepare parameters
        if prec_wrt_L:
            for component in range(nr_components):
                parameter = 'prec_'+str(component)
                parameters_transformed_back[parameter] = list(np.array(parameters_transformed_back[parameter]) * L)

        covMatdiag = np.zeros((nr_components, 400))
        log_det_precMat = np.zeros(nr_components)
        for component in range(nr_components):
            parameter = 'prec_'+str(component)
            covMatdiag[component] = 1.0 / np.array(parameters_transformed_back[parameter])
            log_det_precMat[component] =  np.sum(np.log(parameters_transformed_back[parameter]))

        lik_protein = LikelihoodProtein(braw, Nij, qij)

        lik_protein.set_pairs(
            protein_data['residue_i'],
            protein_data['residue_j'],
            protein_data['contact']
        )

        lik_protein.set_parameters_parallel(
            parameters_transformed_back, covMatdiag, log_det_precMat, nr_components, fixed_parameters, prec_wrt_L)

        #compute f and gradients
        lik_protein.compute_f_df(compute_gradients=compute_gradients)
        f_protein = lik_protein.get_f()
        grad_protein = lik_protein.get_gradients() #={} if no gradients are caluclated

        return f_protein, grad_protein

    # def f_df_protein(self, protein, parameters_structured):
    #     print(protein)
    #
    #     braw = raw.parse_msgpack(self.training_data[protein]['braw_file_path'])
    #     Nij, qij = io.read_qij(self.training_data[protein]['qijabfilename'], braw.ncol)
    #
    #     lik_protein = LikelihoodProtein(braw, Nij, qij)
    #
    #     lik_protein.set_pairs(
    #         self.training_data[protein]['residue_i'],
    #         self.training_data[protein]['residue_j'],
    #         self.training_data[protein]['contact']
    #     )
    #
    #     #create new Parameter instance
    #     parameters = Parameters(self.parameters.parameter_dir)
    #     parameters.set_sigma(self.parameters.sigma, self.parameters.prec_wrt_L)
    #     parameters.set_nr_components(self.parameters.nr_components)
    #     parameters.set_fixed_parameters(self.parameters.fixed_parameters)
    #     parameters.set_parameters_structured(parameters_structured)
    #     lik_protein.set_parameters(parameters)
    #
    #     lik_protein.compute_f_df()
    #
    #     f_protein = lik_protein.get_f()
    #     grad_protein = lik_protein.get_gradients()
    #
    #     return f_protein, grad_protein

    def f_df_py(self, parameters_linear, plot=True, save=True, compute_gradients=True):

        self.parameters.parameters_linear = np.array(parameters_linear)
        self.parameters.parameters_structured = self.parameters.linear_to_structured(parameters_linear)

        #parallel version
        pool = mp.ProcessingPool(self.threads_per_protein)

        #prepare parameters: transform back
        parameters_transformed_back = self.parameters.transform_parameters(weights=True, mean=False, prec=True, back=True)


        f_df_py_protein_parallel = partial(
            self.f_df_py_protein_parallel,
            parameters_transformed_back=parameters_transformed_back,
            prec_wrt_L=self.parameters.prec_wrt_L,
            nr_components=self.parameters.nr_components,
            fixed_parameters=self.parameters.fixed_parameters,
            status=save, #turns on/off printing of protein names,
            compute_gradients=compute_gradients
        )


        #map
        t = time.time()
        results = pool.map(f_df_py_protein_parallel, self.training_data.iteritems())
        f_proteins, grad_proteins = zip(*results)
        timestamp = time.time() - t

        #reduce
        f = np.sum(f_proteins)
        g = {}
        if compute_gradients:
            for parameter in self.parameters.parameters_structured.keys():
                if parameter not in self.parameters.fixed_parameters:
                    g[parameter] = np.sum([grad[parameter] for grad in grad_proteins], axis = 0)



        #sequential version
        #
        #f = 0
        #g = {}
        #
        #initialise gradient dict with zeros
        # for parameter in self.parameters.parameters_structured.keys():
        #     if parameter not in self.parameters.fixed_parameters:
        #         g[parameter] = [0] * len(self.parameters.parameters_structured[parameter])
        #
        #
        # t = time.time()
        # for protein in self.training_data.keys():
        #     f_protein, grad_protein = self.f_df_protein(protein, parameters_structured.copy())
        #     f += f_protein
        #     for parameter, val in grad_protein.iteritems():
        #         g[parameter] += val
        # timestamp = time.time() - t


        # regularization
        reg, reg_gradients = self.compute_regularization(self.parameters.parameters_structured)
        f -= reg
        if compute_gradients:
            for key in reg_gradients.keys():
                g[key] =  list(np.array(g[key]) - np.array(reg_gradients[key]))

        #handle isotrope case: gradients of all 400 dimenstions need to be summed
        if compute_gradients and self.parameters.sigma == 'isotrope':
            for key in g.keys():
                if 'prec' in key:
                    g[key] = [np.sum(g[key])] * 400


        #use mean of function value and gradients
        nr_residue_pairs = self.dataset.nr_pairs_contact + self.dataset.nr_pairs_noncontact
        f /= nr_residue_pairs
        if compute_gradients:
            for key in g.keys():
                g[key] = list(np.array(g[key]) / nr_residue_pairs)


        # save parameters and settings and update log_file
        if save:

            # print log
            self.print_status(f, reg, g, timestamp)

            log_df = self.update_optimization_logfile(f, g, timestamp)

            self.parameters.write_parameters(parameters_transformed_back)
            self.write_settings()

            # Plot the evaluation plot
            if plot:
                opt_plot.plot_evaluation(
                    parameters_transformed_back,
                    log_df,
                    self.get_settings(),
                    self.evaluation_set_kde,
                    self.plot_name
                )

        return(f, self.parameters.structured_to_linear(g))

    def f_df(self, parameters_linear, plot=True, save=True):

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


        #use mean of function value and gradients
        nr_residue_pairs = self.dataset.nr_pairs_contact + self.dataset.nr_pairs_noncontact
        f /= nr_residue_pairs
        for key in gradients.keys():
            gradients[key] = list(np.array(gradients[key]) / nr_residue_pairs)


        # save parameters and settings and update log_file
        if save:

            # print log
            self.print_status(f, reg, gradients, timestamp)

            log_df = self.update_optimization_logfile(f, gradients, timestamp)

            parameters_transformed_back = self.parameters.transform_parameters(weights=True, mean=False, prec=True, back=True)
            self.parameters.write_parameters(parameters_transformed_back)
            self.write_settings()

            # Plot the evaluation plot
            if plot:
                opt_plot.plot_evaluation(
                    parameters_transformed_back,
                    log_df,
                    self.get_settings(),
                    self.evaluation_set_kde,
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

        #compute mean of function value for cross-val set
        nr_residue_pairs = self.dataset.nr_pairs_contact_cross_val + self.dataset.nr_pairs_noncontact_cross_val
        f_crossval /= nr_residue_pairs


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

