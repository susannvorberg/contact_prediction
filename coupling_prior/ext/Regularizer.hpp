//============================================================================
// Name        : Regularizer.hpp
// Author      : susi
// Version     :
// Copyright   : Your copyright notice
// Description : header file for regularizer of ll
//============================================================================

#ifndef REGULARIZER
#define REGULARIZER

#include "Parameters.hpp"
#include <armadillo>
#include <map>

/********************************************************************************************************************************
 *
 * Class header
 *
 *********************************************************************************************************************************/

class Regularizer{
	public:
		//constructor
		Regularizer(std::map<std::string, std::vector<double> > parameterMap,
		            double regularization_parameter_mu_,
                    double regularization_parameter_diagonal_PrecMat_,
		            int debug_
		            );
		//default constructor
		Regularizer();


		/*
		 * Set the regularization parameters prior to calculation of the regularization term
		 */
		void set_regularization_parameters(double regularization_parameter_mu_, double regularization_parameter_diagonal_PrecMat_);


		/*
		 * Print regularization coefficients to stdout
		 */
		void print_regularization_parameters();


		/**********************************************************************
		 * Compute regularizers
		 **********************************************************************/
		/*
		 * Gaussian prior on Mean parameter
		 *
		 * --> penalizes deviation from Zero centered Gaussian
		 *
		 * N(mu | 0, "regularization_parameter_mu")
		*/
		double regularizer_mu();


		/*
		 * Gaussian prior on diagonal of precision matrix Lambda_k
		 *
		 * --> penalizing values that are too precise == very sharp Gaussians
		 *
		 * N(diag(Lambda)| 0, "regularization_diagonal_covMat")
		 */
		double regularizer_diagPrecMat();

		/**********************************************************************
		 * Gradients of regularizers
		 **********************************************************************/

		/*
		 * This function calculates the gradient of the regularizer of mu
		 * for a component c
		 */
		arma::vec gradient_mu_comp_reg(int k);

		/*
		 * Calculate the gradient of the regularizer of the diagonal precision matrix
		 * for a component c
		 */
		//arma::mat gradient_choleskyL_comp_reg(int c);
        arma::mat gradient_diagprecMat_comp_reg(int k);


	private:
		int debug;
		int nr_components;
	    Parameters parameters;
		double regularization_parameter_mu;
		double regularization_parameter_diagonal_PrecMat;

};

#endif //REGULARIZER
