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

        str = "Hyperparameters for Coupling Prior will be inferred using the following settings:\n"

        str += self.likelihood.dataset.__repr__()
        str += self.likelihood.__repr__()

        str += "\nOptimization specific settings: \n"
        str += "maxiter: \t\t\t\t\t\t {0} \n".format(self.maxiter)
        str += "method: \t\t\t\t\t\t {0} \n".format(self.method)


        return(str)

    def minimize(self):

        self.likelihood.debug_mode = self.debug_mode

        res = scipy.optimize.minimize(
            self.likelihood.f_df,
            self.likelihood.get_parameters_linear(),
            method=self.method,
            jac=True,
            options={
                'maxiter': self.maxiter
                , 'disp': True
            }
        )
        print(res)