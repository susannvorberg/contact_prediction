#!/usr/bin/env python

import scipy.optimize


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

        str += self.likelihood.dataset.__repr__()
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

    def minimize(self):

        self.likelihood.debug_mode = self.debug_mode

        res = scipy.optimize.minimize(
            self.likelihood.f_df,
            self.likelihood.get_parameters_linear(),
            method=self.method,
            jac=True,
            options={
                'maxiter': self.maxiter
                , 'disp': False
            }
        )
        print(res)

        parameters_struct = self.likelihood.linear_to_structured(res.x)

