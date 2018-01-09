#!/usr/bin/env python


import numpy as np


class Regularization():
    """
    Specify the Regularization of Likelihood Parameters for Copupling Prior
    """

    def __init__(self, parameters, reg_coeff_mu, reg_coeff_diag_prec):

        self.fixed_parameters = parameters.fixed_parameters
        self.parameters_mu = {}
        self.parameters_diag_prec = {}
        self.set_parameters(parameters.parameters_structured)

        #regularization coefficient
        self.reg_coeff_mu          = 0
        self.reg_coeff_diag_prec   = 0
        self.set_reg_coefficients(reg_coeff_mu, reg_coeff_diag_prec)

    def __repr__(self):

        str = "\nRegularization settings: \n"
        for param in ["self.reg_coeff_mu",
                      "self.reg_coeff_diag_prec"]:
            str += "{0:{1}} {2}\n".format(param.split(".")[1], '36', eval(param))

        return str

    def set_parameters(self, parameters_structured):
        for parameter_name, parameter in parameters_structured.iteritems():
            if parameter_name not in self.fixed_parameters:

                if 'mu' in parameter_name:
                    self.parameters_mu[parameter_name] = parameter

                if 'prec' in parameter_name:
                    #transform
                    self.parameters_diag_prec[parameter_name] = np.exp(parameter)

    def set_reg_coefficients(self, reg_coeff_mu, reg_coeff_diag_prec):
        self.reg_coeff_mu                 = reg_coeff_mu
        self.reg_coeff_diag_prec   = reg_coeff_diag_prec

    def get_settings(self):
        settings={}
        settings['reg_coeff_mu'] = self.reg_coeff_mu
        settings['reg_coeff_diag_prec'] = self.reg_coeff_diag_prec

        return settings

    def get_regularization(self):
        return self._regularizer_diag_prec() + self._regularizer_mu()

    def get_regularization_gradients(self):
        gradients = {}

        if self.reg_coeff_mu != 0:
            for parameter in self.parameters_mu.keys():
                gradients[parameter] = -np.array(self.parameters_mu[parameter]) / (self.reg_coeff_mu*self.reg_coeff_mu)
        else:
            for parameter in self.parameters_mu.keys():
                gradients[parameter] = [0] * len(self.parameters_mu[parameter])


        if self.reg_coeff_diag_prec != 0:
            for parameter in self.parameters_diag_prec.keys():
                gradients[parameter] = -np.array(self.parameters_diag_prec[parameter])  / (self.reg_coeff_diag_prec * self.reg_coeff_diag_prec)
                # add derivative of exponential transformation
                gradients[parameter] *= np.array(self.parameters_diag_prec[parameter])

                #in case of regularizing the covariance matrix,
                #the gradient looks like (inkl derivative of exp):
                #gradients[parameter] = 1/(self.reg_coeff_diag_prec * self.reg_coeff_diag_prec) * 1/(np.array(self.parameters_diag_prec[parameter]) * np.array(self.parameters_diag_prec[parameter]))

        else:
            for parameter in self.parameters_diag_prec.keys():
                gradients[parameter] = [0] * len(self.parameters_diag_prec[parameter])

        return gradients

    def _regularizer_mu(self):
        reg = 0
        if self.reg_coeff_mu != 0:

            for parameter in self.parameters_mu.values():
                reg += np.sum(np.array(parameter) * np.array(parameter))

            #defines regularization strength
            reg *= -1.0/(2*self.reg_coeff_mu*self.reg_coeff_mu)

        return reg

    def _regularizer_diag_prec(self):
        reg = 0
        if self.reg_coeff_diag_prec != 0:

            for parameter in self.parameters_diag_prec.values():
                reg += np.sum(np.array(parameter) * np.array(parameter))

                #in case of regularizing the covariance matrix,
                #the regularizer looks as follows:
                #reg += np.sum(np.array(1/parameter) * np.array(1/parameter))

            #defines regularization strength
            reg *= -1.0/(2*self.reg_coeff_diag_prec*self.reg_coeff_diag_prec)

        return reg



