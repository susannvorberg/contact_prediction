import numpy as np

class LikelihoodProtein():
    """
    Specify the Likelihood and Gradients

    """

    def __init__(self, braw, Nij, qij ):

        self.braw   = braw
        self.wij = self.braw.x_pair[:,:,:20, :20].reshape((braw.ncol, braw.ncol, 400))
        self.Nij    = Nij
        self.qij    = qij

        self.L      = braw.ncol

        self.lambda_w =  braw.meta['workflow'][0]['regularization']['lambda_pair']
        self.lambda_w_mat = np.diag([self.lambda_w] * 400)

        self.residues_i = None
        self.residues_j = None
        self.contacts = None

        self.parameters = None
        self.covMatdiag = None
        self.log_det_precMat = None

        self.f = 0

        #intialize gradients
        self.gradients = {}


    def set_pairs(self, i,j,contact):
        self.residues_i = i
        self.residues_j = j
        self.contacts = contact

    def set_parameters(self, parameters):

        self.parameters = parameters
        self.covMatdiag = np.zeros((self.parameters.nr_components, 400))
        self.log_det_precMat = np.zeros(self.parameters.nr_components)

        for parameter in self.parameters.parameters_structured.keys():
            if parameter not in self.parameters.fixed_parameters:
                self.gradients[parameter] = [0] * len(self.parameters.parameters_structured[parameter])

                if 'weight' in parameter:
                    #transform
                    self.parameters.parameters_structured[parameter] = np.exp(self.parameters.parameters_structured[parameter])

                if 'prec' in parameter:
                    component = int(parameter.split("_")[1])

                    #transform
                    self.parameters.parameters_structured[parameter] = np.exp(self.parameters.parameters_structured[parameter])

                    if self.parameters.prec_wrt_L:
                        self.parameters.parameters_structured[parameter] *= self.L

                    self.covMatdiag[component] = 1.0 / np.array(self.parameters.parameters_structured[parameter])
                    self.log_det_precMat[component] =  np.sum(np.log(self.parameters.parameters_structured[parameter]))

    def get_f(self):
        return self.f

    def get_gradients(self):
        return self.gradients

    def log_density_gaussian_ratio(self, log_det_lambdak, mu_k, lambda_k, log_det_lambdaijk, mu_ij_k, lambda_ij_k):

        gaussian_1 = log_det_lambdak   - np.matmul(mu_k,    np.matmul(lambda_k,    mu_k))
        gaussian_2 = log_det_lambdaijk - np.matmul(mu_ij_k, np.matmul(lambda_ij_k, mu_ij_k))

        gaussian_ratio_log_density = 0.5 * (gaussian_1 - gaussian_2)

        return gaussian_ratio_log_density

    def compute_gradients_protein(self, contact, responsibilities, mu_ij_k, lambda_ij_k_inv):

        for parameter in self.parameters.parameters_structured.keys():
            if parameter not in self.parameters.fixed_parameters:

                component = int(parameter.split("_")[-1])

                grad = [0]
                if "weight" in parameter:
                    grad = self.compute_gradient_weight(contact, component, responsibilities)
                elif "prec" in parameter:
                    grad = self.compute_gradient_precMat(component, mu_ij_k, lambda_ij_k_inv, responsibilities)
                elif "mu" in parameter:
                    grad = self.compute_gradient_mu(component, mu_ij_k, responsibilities)

                self.gradients[parameter] += grad

    def compute_gradient_weight(self, contact, component, responsibilities ):
        """

        :param contact:
        :param component:
        :param responsibilities:
        :return: gradient for LOG LIKELIHOOD
        """
        if contact:
            weight_k = self.parameters.parameters_structured['weight_contact_'+str(component)]
        else:
            weight_k = self.parameters.parameters_structured['weight_bg_'+str(component)]

        grad_weight = responsibilities[component] - weight_k

        #gradient of NEG log likelihood
        return -grad_weight

    def compute_gradient_mu(self, component, mu_ij_k, responsibilities):

        mu_k = self.parameters.parameters_structured['mu_'+str(component)]
        precMatdiag = self.parameters.parameters_structured['prec_'+ str(component)]

        diff = mu_ij_k[component] - mu_k

        #precMatdiag *  diff => element-wise because precMat is always diagonal
        grad_mu = responsibilities[component] * (precMatdiag *  diff)

        #gradient of NEG log likelihood
        return -grad_mu

    def compute_gradient_precMat(self, component, mu_ij_k, lambda_ij_k_inv, responsibilities):
        mu_k = self.parameters.parameters_structured['mu_'+str(component)]
        precMatdiag = np.array(self.parameters.parameters_structured['prec_'+ str(component)])
        covMatdiag = np.diag(self.covMatdiag[component])

        diff = mu_ij_k[component] - mu_k
        outer_product_diff = np.outer(diff, diff)

        grad_precMat = 0.5 * responsibilities[component] * (covMatdiag - lambda_ij_k_inv[component] - outer_product_diff )

        #derivative of exponential transformation
        grad_precMat[np.diag_indices_from(grad_precMat)] *= precMatdiag


        if self.parameters.sigma == 'isotrope':
            if self.parameters.prec_wrt_L:
                grad_precMat *= self.L
            grad = np.full(400, -np.sum(np.diag(grad_precMat)))
        elif self.parameters.sigma == 'diagonal':
            grad = np.diag(-grad_precMat)
        else:
            grad = -grad_precMat

        #gradient of NEG log likelihood
        return grad

    def compute_f_df(self, compute_gradients=True):

        if self.parameters is None:
            print("You first need to set parameters!")
            return

        #initialise f and gradients in case compute_f_df is called multiple times
        self.f = 0
        for parameter in self.parameters.parameters_structured.keys():
            if parameter not in self.parameters.fixed_parameters:
                self.gradients[parameter] = [0] * len(self.parameters.parameters_structured[parameter])

        #iterate over all residue pairs for this protein
        for i,j,contact in zip(self.residues_i, self.residues_j, self.contacts):


            diag_qij 		= np.diag(self.qij[i,j])                    # dim 400, 400 only diagonal
            qij_prod        = np.outer(self.qij[i,j], self.qij[i,j])    # dim 400, 400 = matrix multiplication u * vT == outer product of two vectors

            Hij = self.Nij[i,j] * (diag_qij - qij_prod)  + self.lambda_w_mat #eq 37       # dim 400, 400
            Hij_wij_prod = np.matmul(Hij, self.wij[i,j])                     #[400x400] * [400] => [400]


            # print   "i: " + str(i) + " j: " + str(j) + " contact:  " + str(contact)
            # print   "vqij(21): " + str(self.qij[i,j,21]) + " vqij(22): " + str(self.qij[i,j,22])
            # print   "N_ij: " + str(self.Nij[i,j])
            # print   "w_ij(0): " + str(self.wij[i,j,0] ) +  " w_ij(1): " + str(self.wij[i,j,1]) +  " w_ij(2): " + str(self.wij[i,j,2])
            # print   "qij_prod(0,1): " + str(qij_prod[0,1])  +  " qij_prod(0,2): " + str(qij_prod[0,2] )+  " qij_prod(2,0): " + str(qij_prod[2,0])
            # print   "H_ij(0,1): " + str(Hij[0,1])  +  " H_ij(0,2): " + str(Hij[0,2]) +  " H_ij(2,0): " + str(Hij[2,0])
            # print   "Hij_wij_prod(0): " + str(Hij_wij_prod[0])  +  " Hij_wij_prod(1): " + str(Hij_wij_prod[1]) +  " Hij_wij_prod(2): " + str(Hij_wij_prod[2])

            lambda_ij_k = np.zeros((self.parameters.nr_components, 400, 400))
            lambda_ij_k_inv = np.zeros((self.parameters.nr_components, 400, 400))
            mu_ij_k = np.zeros((self.parameters.nr_components, 400))
            log_det_lambda_ij_k = np.zeros(self.parameters.nr_components)
            log_density = np.zeros(self.parameters.nr_components)

            for component in range(self.parameters.nr_components):

                #define component specific parameters
                precMat = np.diag(self.parameters.parameters_structured['prec_'+ str(component)])
                precMatdiag = np.array(self.parameters.parameters_structured['prec_'+ str(component)])
                mu_k = self.parameters.parameters_structured['mu_'+str(component)]
                if contact:
                    weight_k = self.parameters.parameters_structured['weight_contact_'+str(component)]
                else:
                    weight_k = self.parameters.parameters_structured['weight_bg_'+str(component)]


                #---------------- simplify computation in case that lambda_k is diagonal matrix
                A = self.Nij[i,j] * self.qij[i,j] + precMatdiag  #represents diagonal of 400d matrix
                A_inv = 1.0 / A                             #represents inverse of diagonal 400d matrix
                log_det_A	= np.sum(np.log(A))
                Ainv_qij_product  = A_inv * self.qij[i,j]        #element-wise multiplication (diagMat * vector => vector)
                triple_product 	    = np.sum(self.qij[i,j] * Ainv_qij_product) #element-wise multiplication (vector^T * diagMat * vector = scalar)
                A_qij_outer_product = np.outer(Ainv_qij_product, Ainv_qij_product)

                #compute lambda_ij_k
                lambda_ij_k[component]      = Hij - self.lambda_w_mat + precMat

                lambda_ij_k_inv[component]  = np.diag(A_inv) + A_qij_outer_product / (1.0/self.Nij[i,j] - triple_product)


                #compute mu_ij_k
                mu_ij_k[component]   = np.matmul(lambda_ij_k_inv[component], ( Hij_wij_prod + (precMatdiag * mu_k)))


                #save log determinant of lambda_ij_k, see page 16 Saikats theory
                log_det_lambda_ij_k[component] = np.log(1 - self.Nij[i,j] * triple_product) + log_det_A


                #ratio of two gaussian log densities
                gaussian_ratio_logdensity = self.log_density_gaussian_ratio(
                    self.log_det_precMat[component], mu_k, precMat,
                    log_det_lambda_ij_k[component], mu_ij_k[component], lambda_ij_k[component]
                )


                log_density[component] = np.log(weight_k) + gaussian_ratio_logdensity


                # if component ==0:
                #     print  "triple_product: " + str(triple_product)
                #     print  "log_det_A: " + str(log_det_A)
                #     print  "lambda_ij_k_mat(0,1): " + str(lambda_ij_k[component,0,1])  +  " lambda_ij_k_mat(0,2): " + str(lambda_ij_k[component,0,2]) +  " lambda_ij_k_mat(2,0): " + str(lambda_ij_k[component,2,0])
                #     print  "Gaussian log density : " + str(gaussian_ratio_logdensity)
                #     print  "mu_ij_k_vec(0): " + str(mu_ij_k[component,0]) +  " mu_ij_k_vec(1): " + str(mu_ij_k[component,1]) +  " mu_ij_k_vec(2): " + str(mu_ij_k[component,2])




            #Johannes suggestion how to precisely compute the responsibilities
            a_max = np.max(log_density)
            responsibilities = np.exp(log_density - a_max)
            sum_resp = np.sum(responsibilities)
            responsibilities /= sum_resp



            #save NEG log likelihood of current pair
            f = np.log(sum_resp) + a_max
            if not np.isnan(f) and not np.isinf(f):
                self.f -= f

            # print   "resp: " + str(responsibilities[0]) + " f: " + str(f)

            if compute_gradients:
                self.compute_gradients_protein(contact, responsibilities, mu_ij_k, lambda_ij_k_inv)



