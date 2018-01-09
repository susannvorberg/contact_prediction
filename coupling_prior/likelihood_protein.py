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


        if 'regularization' in braw.meta['workflow'][0]:
            self.lambda_w =  braw.meta['workflow'][0]['regularization']['lambda_pair']
        else:
            self.lambda_w = braw.meta['workflow'][0]['parameters']['regularization']['lambda_pair']


        self.lambda_w_mat = np.diag([self.lambda_w] * 400)

        self.residues_i = None
        self.residues_j = None
        self.contacts = None

        self.f_pairwise = []
        #intialize gradients
        self.gradients = {}


        self.nr_components = None
        self.parameters_structured = None
        self.fixed_parameters = None
        self.prec_wrt_L=False
        self.covMatdiag = None
        self.log_det_precMat = None

    def set_pairs(self, i,j,contact):
        self.residues_i = i
        self.residues_j = j
        self.contacts = contact

    def set_parameters_parallel(self, parameters, covMatdiag, log_det_precMat, nr_components, fixed_parameters, prec_wrt_L):

        self.nr_components = nr_components
        self.parameters_structured = parameters
        self.fixed_parameters = fixed_parameters
        self.prec_wrt_L=prec_wrt_L
        self.covMatdiag = covMatdiag
        self.log_det_precMat = log_det_precMat

        #initialise gradient dict with zeros - only for parameters that will be optimized!
        for parameter in self.parameters_structured.keys():
            if parameter not in self.fixed_parameters:
                self.gradients[parameter] = [0] * len(self.parameters_structured[parameter])


    def set_parameters(self, parameters):


        self.parameters = parameters


        sum_weight_contact = 0
        sum_weight_bg = 0
        for component in range(self.parameters.nr_components):

            parameter = 'weight_contact_'+str(component)
            if parameter not in self.parameters.fixed_parameters:
                self.gradients[parameter] = [0] * len(self.parameters.parameters_structured[parameter])
            self.parameters.parameters_structured[parameter] = np.exp(self.parameters.parameters_structured[parameter])
            sum_weight_contact += self.parameters.parameters_structured[parameter]


            parameter = 'weight_bg_'+str(component)
            if parameter not in self.parameters.fixed_parameters:
                self.gradients[parameter] = [0] * len(self.parameters.parameters_structured[parameter])
            self.parameters.parameters_structured[parameter] = np.exp(self.parameters.parameters_structured[parameter])
            sum_weight_bg += self.parameters.parameters_structured[parameter]


            parameter = 'mu_'+str(component)
            if parameter not in self.parameters.fixed_parameters:
                self.gradients[parameter] = [0] * len(self.parameters.parameters_structured[parameter])


            parameter = 'prec_'+str(component)
            if parameter not in self.parameters.fixed_parameters:
                self.gradients[parameter] = [0] * len(self.parameters.parameters_structured[parameter])
            self.parameters.parameters_structured[parameter] = np.exp(self.parameters.parameters_structured[parameter])

            if self.parameters.prec_wrt_L:
                self.parameters.parameters_structured[parameter] *= self.L

            self.covMatdiag[component] = 1.0 / np.array(self.parameters.parameters_structured[parameter])
            self.log_det_precMat[component] =  np.sum(np.log(self.parameters.parameters_structured[parameter]))


        for component in range(self.parameters.nr_components):
            parameter = 'weight_contact_'+str(component)
            self.parameters.parameters_structured[parameter] /= sum_weight_contact

            parameter = 'weight_bg_'+str(component)
            self.parameters.parameters_structured[parameter] /= sum_weight_bg

    def get_f(self):
        return np.sum(self.f_pairwise)

    def get_f_pairwise(self):
        return self.f_pairwise

    def get_gradients(self):
        return self.gradients

    @staticmethod
    def log_density_gaussian_ratio(log_det_lambdak, mu_k, lambda_k, log_det_lambdaijk, mu_ij_k, lambda_ij_k):

        gaussian_1 = log_det_lambdak   - np.matmul(mu_k,    np.matmul(lambda_k,    mu_k))
        gaussian_2 = log_det_lambdaijk - np.matmul(mu_ij_k, np.matmul(lambda_ij_k, mu_ij_k))

        gaussian_ratio_log_density = 0.5 * (gaussian_1 - gaussian_2)

        return gaussian_ratio_log_density

    def compute_gradients_protein(self, contact, responsibilities, mu_ij_k, lambda_ij_k_inv):


        for parameter in self.gradients.keys():
                component = int(parameter.split("_")[-1])

                if "weight_bg" in parameter and not contact:
                    self.gradients[parameter] += self.compute_gradient_weight(contact, component, responsibilities)
                elif "weight_contact" in parameter and contact:
                    self.gradients[parameter] += self.compute_gradient_weight(contact, component, responsibilities)
                elif "prec" in parameter:
                    self.gradients[parameter] += self.compute_gradient_precMat(component, mu_ij_k, lambda_ij_k_inv, responsibilities)
                elif "mu" in parameter:
                    self.gradients[parameter] += self.compute_gradient_mu(component, mu_ij_k, responsibilities)

    def compute_gradient_weight(self, contact, component, responsibilities ):
        """

        :param contact:
        :param component:
        :param responsibilities:
        :return: gradient for LOG LIKELIHOOD
        """

        if contact:
            weight_k = self.parameters_structured['weight_contact_'+str(component)]
        else:
            weight_k = self.parameters_structured['weight_bg_'+str(component)]

        grad_weight = responsibilities[component] - weight_k

        #gradient of NEG log likelihood
        return -grad_weight

    def compute_gradient_mu(self, component, mu_ij_k, responsibilities):

        mu_k = self.parameters_structured['mu_'+str(component)]
        precMatdiag = self.parameters_structured['prec_'+ str(component)]

        diff = mu_ij_k[component] - mu_k

        #precMatdiag *  diff => element-wise because precMat is always diagonal
        grad_mu = responsibilities[component] * (precMatdiag *  diff)

        #gradient of NEG log likelihood
        return -grad_mu

    def compute_gradient_precMat(self, component, mu_ij_k, lambda_ij_k_inv, responsibilities):
        mu_k = self.parameters_structured['mu_'+str(component)]
        precMatdiag = np.array(self.parameters_structured['prec_'+ str(component)])
        covMatdiag = np.diag(self.covMatdiag[component])

        diff = mu_ij_k[component] - mu_k
        outer_product_diff = np.outer(diff, diff)

        grad_precMat = 0.5 * responsibilities[component] * (covMatdiag - lambda_ij_k_inv[component] - outer_product_diff )

        #derivative of exponential transformation (includes derivative of L if self.prec_wrt_L=True)
        grad_precMat[np.diag_indices_from(grad_precMat)] *= precMatdiag

        #derivative of L
        #if self.prec_wrt_L:
        #    grad_precMat *= self.L

        #diagonal of precision matrix
        grad = np.diag(-grad_precMat)

        #gradient of NEG log likelihood
        return grad

    def compute_f_df(self, compute_gradients=True):

        if self.parameters_structured is None:
            print("You first need to set parameters with 'set_parameters()'!")
            return

        if self.residues_i is None:
            print("You first need to set data with 'set_pairs()'!")
            return

        #initialise f and gradients in case compute_f_df is called multiple times
        self.f_pairwise = [0] * len(self.residues_i)
        # for parameter in self.parameters_structured.keys():
        #     if parameter not in self.parameters.fixed_parameters:
        #         self.gradients[parameter] = [0] * len(self.parameters.parameters_structured[parameter])


        #iterate over all residue pairs for this protein
        for nr in range(len(self.residues_i)):

            i = self.residues_i[nr]
            j = self.residues_j[nr]
            contact = self.contacts[nr]


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

            lambda_ij_k = np.zeros((self.nr_components, 400, 400))
            lambda_ij_k_inv = np.zeros((self.nr_components, 400, 400))
            mu_ij_k = np.zeros((self.nr_components, 400))
            log_det_lambda_ij_k = np.zeros(self.nr_components)
            log_density = np.zeros(self.nr_components)

            for component in range(self.nr_components):


                #define component specific parameters
                precMat = np.diag(self.parameters_structured['prec_'+ str(component)])
                precMatdiag = np.array(self.parameters_structured['prec_'+ str(component)])
                mu_k = self.parameters_structured['mu_'+str(component)]
                if contact:
                    weight_k = self.parameters_structured['weight_contact_'+str(component)]
                else:
                    weight_k = self.parameters_structured['weight_bg_'+str(component)]


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

                # if i == 0 and j == 12:
                #     print("\ni={0} j={1} contact={2}".format(i,j,contact))
                #     print  "triple_product: " + str(triple_product)
                #     print  "log_det_A: " + str(log_det_A)
                #     print  "lambda_ij_k_mat(0,1): " + str(lambda_ij_k[component,0,1])  +  " lambda_ij_k_mat(0,2): " + str(lambda_ij_k[component,0,2]) +  " lambda_ij_k_mat(2,0): " + str(lambda_ij_k[component,2,0])
                #     print  "mu_ij_k_vec(0): " + str(mu_ij_k[component,0]) +  " mu_ij_k_vec(1): " + str(mu_ij_k[component,1]) +  " mu_ij_k_vec(2): " + str(mu_ij_k[component,2])
                #     print  "Gaussian log density : " + str(gaussian_ratio_logdensity)
                #     print  "log_density[component]: " + str(log_density[component])
                #     print  "np.log(weight_k) : " + str(np.log(weight_k)) + " weight_k: " + str(weight_k)




            #Johannes suggestion how to precisely compute the responsibilities
            a_max = np.max(log_density)
            responsibilities = np.exp(log_density - a_max)
            sum_resp = np.sum(responsibilities)
            responsibilities /= sum_resp

            #save NEG log likelihood of current pair
            f = np.log(sum_resp) + a_max
            if not np.isnan(f) and not np.isinf(f):
                self.f_pairwise[nr] = -f

            # if i == 0 and j == 12:
            #     print  "a_max: " + str(a_max)
            #     print  "sum_resp: " + str(sum_resp)
            #     print("i={0} j={1} f={2}".format(i,j,f))


            # print   "resp: " + str(responsibilities[0]) + " f: " + str(f)

            if compute_gradients:
                self.compute_gradients_protein(contact, responsibilities, mu_ij_k, lambda_ij_k_inv)



