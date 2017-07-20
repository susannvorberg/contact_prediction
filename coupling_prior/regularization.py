#!/usr/bin/env python


import numpy as np


class Regularization():
    """
    Specify the Regularization of Likelihood Parameters for Copupling Prior
    """

    def __init__(self, parameters, reg_coeff_mu, reg_coeff_diag_prec):

        self.parameters_mu = {}
        self.parameters_diag_prec = {}
        self.set_parameters(parameters)

        #regularization coefficient
        self.reg_coeff_mu          = 0
        self.reg_coeff_diag_prec   = 0
        self.set_reg_coefficients(reg_coeff_mu, reg_coeff_diag_prec)

    def set_parameters(self, parameters):
        for parameter_name, parameter in parameters.parameters_structured.keys():
            if parameter_name not in parameters.fixed_parameters:

                if 'mu' in parameter_name:
                    self.parameters_mu[parameter_name] = parameter

                if 'prec' in parameter:
                    #transform
                    self.parameters_diag_prec[parameter_name] = np.exp(parameter)

    def set_reg_coefficients(self, reg_coeff_mu, reg_coeff_diag_prec):
        self.reg_coeff_mu                 = reg_coeff_mu
        self.reg_coeff_diag_prec   = reg_coeff_diag_prec

    def get_regularization(self):
        return self._regularizer_diag_prec() + self._regularizer_mu()

    def get_regularization_gradients(self):
        gradients = {}

        for parameter in self.parameters_mu.keys():
            gradients[parameter] = - self.parameters_mu[parameter] / (self.reg_coeff_mu*self.reg_coeff_mu)

        for parameter in self.parameters_diag_prec.keys():
            gradients[parameter] = -self.parameters_diag_prec[parameter] / (self.reg_coeff_diag_prec * self.reg_coeff_diag_prec)
            # add derivative of exponential transformation
            gradients[parameter] *= self.parameters_diag_prec[parameter]

        return gradients

    def _regularizer_mu(self):
        reg = 0
        if self.reg_coeff_mu != 0:

            for parameter in self.parameters_mu.values():
                reg += np.linalg.norm(parameter, ord=2)

            #defines regularization strength
            reg *= -1/(2*self.reg_coeff_mu*self.reg_coeff_mu)

        return reg

    def _regularizer_diag_prec(self):
        reg = 0
        if self.reg_coeff_diag_prec != 0:
            for parameter in self.parameters_diag_prec:
                reg += np.linalg.norm(parameter, ord=2)

            #defines regularization strength
            reg *= -1/(2*self.reg_coeff_diag_prec*self.reg_coeff_diag_prec)

        return reg



