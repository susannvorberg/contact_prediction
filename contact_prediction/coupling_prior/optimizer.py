#!/usr/bin/env python

import scipy.optimize
import json

class Optimizer():
    """
    Optimize Coupling Prior with specified Optimizer
    """
    def __init__(self, likelihood):


        self.likelihood = likelihood

        #optimization settings
        self.maxiter                = 1000
        self.method                 = 'L-BFGS-B' #'CG'
        self.debug_mode             = 0

    def __repr__(self):

        str = "\nHyperparameters for Coupling Prior will be inferred using the following settings:\n"

        str += self.likelihood.__repr__()

        str += "\nOptimization specific settings: \n"
        str += "{0:{1}} {2}\n".format("max nr iterations", '36', self.maxiter)
        str += "{0:{1}} {2}\n".format("optimizer", '36', self.method)

        return(str)

    def set_maxiter(self, maxiter):
        self.maxiter = int(maxiter)

    def set_method(self, method):
        if method not in ['L-BFGS-B', 'CG']:
            self.method = 'L-BFGS-B'
            print("Method must be one of ['L-BFGS-B', 'CG']!")
        else:
            self.method = method

    def set_debug_mode(self, debug_mode):
        self.debug_mode = int(debug_mode)

    def write_res(self, settings, res, settings_file):

        settings['result']={}
        settings['result']['f_value'] = res.fun
        settings['result']['nr_iterations'] = res.nit
        settings['result']['nr_evaluations'] = res.nfev
        settings['result']['success'] = res.success
        settings['result']['opt_message'] = res.message

        with open(settings_file, 'w') as fp:
            json.dump(settings, fp)

    def minimize(self):

        self.likelihood.debug_mode = self.debug_mode

        f_df = self.likelihood.f_df
        if self.likelihood.python_parallel:
            f_df= self.likelihood.f_df_py

        res = scipy.optimize.minimize(
            f_df,
            self.likelihood.parameters.get_parameters_linear(),
            method=self.method,
            jac=True,
            options={
                'maxiter': self.maxiter
                , 'disp': False
            }
        )
        print(res)

        self.likelihood.parameters.parameters_linear = res.x
        self.likelihood.parameters.parameters_structured = self.likelihood.parameters.linear_to_structured(res.x)
        parameters_transformed_back = self.likelihood.parameters.transform_parameters(weights=True, mean=False, prec=True, back=True)

        ##### save parameters and settings
        self.likelihood.parameters.write_parameters(parameters_transformed_back)
        self.write_res(self.likelihood.get_settings(), res, self.likelihood.settings_file)

