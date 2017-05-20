//============================================================================
// Name        : Regularizer.cpp
// Author      : Susann Vorberg
// Version     :
// Copyright   : Your copyright notice
// Description : implementation of regularizer for parameters of ll
//============================================================================

#include "Regularizer.hpp"
#include "Parameters.hpp"
#include <armadillo>
#include <map>

/********************************************************************************************************************************
 *
 * Class constructor
 *
 * @param parameterMap_
 * @param debug_mode
 *
 *
 *********************************************************************************************************************************/


Regularizer::Regularizer(   std::map<std::string, std::vector<double> > parameterMap,
                            double regularization_parameter_mu_,
                            double regularization_parameter_diagonal_PrecMat_
                         ):parameters(parameterMap.size()/ 4)
{

    //set up model specifications
    this->regularization_parameter_mu 	= regularization_parameter_mu_;
	this->regularization_parameter_diagonal_PrecMat 	= regularization_parameter_diagonal_PrecMat_;
	this->nr_components = parameters.get_nr_components();


	//set up parameters
    parameters.set_parameters(parameterMap);

}




/********************************************************************************************************************************
 *
 * Default Class constructor
 *
 *********************************************************************************************************************************/

Regularizer::Regularizer()
{
}



/*
 * Print regularization coefficients to stdout
 */
void Regularizer::print_regularization_parameters(){

	std::cout << "\nRegularization coefficient mu: " << this->regularization_parameter_mu  << std::endl;
	std::cout << "Regularization coefficient diagonal Precision matrix:  " << this->regularization_parameter_diagonal_PrecMat  << std::endl;

}




/*
 * Set the regularization parameters prior to calculation of the regularization term
 */
void Regularizer::set_regularization_parameters(double regularization_parameter_mu_, double regularization_parameter_diagonal_PrecMat_) {

	this->regularization_parameter_mu 	= regularization_parameter_mu_;
	this->regularization_parameter_diagonal_PrecMat 	= regularization_parameter_diagonal_PrecMat_;
}



/********************************************************************************************************************************
 *
 * Regualarizers
 *
 *********************************************************************************************************************************/

/*
 * Gaussian prior on Mean parameter
 *
 * --> penalizes deviation from Zero centered Gaussian
 *
 * N(mu | 0, "regularization_parameter_mu")
*/
double Regularizer::regularizer_mu(){

	double reg 	= 0.0;

	for(int k = 1; k < this->nr_components; k++){
		reg += arma::accu(arma::square(this->parameters.get_mean(k)));
	}

	reg *= -1/(2 * this->regularization_parameter_mu * this->regularization_parameter_mu);

	return(reg);
}

/*
 * Gaussian prior on diagonal of covariance matrix Lambda_k^-1
 *
 * --> penalizing values that are too broad == very flat Gaussians
 *
 * N(diag(Lambda^-1)| 0, "regularization_diagonal_covMat")
 */
double Regularizer::regularizer_diagPrecMat(){

	double reg 	= 0.0;
	for(int k = 1; k < this->nr_components; k++){
		arma::mat precMat_k   = this->parameters.get_precMat(k);
		reg += arma::accu(arma::square(precMat_k.diag())); // square for l2-regularization
	}
	reg 	*= -1/(2 * this->regularization_parameter_diagonal_PrecMat * this->regularization_parameter_diagonal_PrecMat);

	return(reg);
}

/********************************************************************************************************************************
 *
 * Gradients of Regularizers
 *
 *********************************************************************************************************************************/


/*
 * This function calculates the gradient of the regularizer of mu
 */
arma::vec Regularizer::gradient_mu_comp_reg(int k){

	arma::vec mean 		= this->parameters.get_mean(k);
	arma::vec grad_mean = mean / (this->regularization_parameter_mu * this->regularization_parameter_mu);

	return (-grad_mean);

}


/*
 * Calculate the gradient of the regularizer of the diagonal precision matrix
 *
 */
arma::mat Regularizer::gradient_diagprecMat_comp_reg(int k){

	arma::mat precMat_k 	= this->parameters.get_precMat(k);

    arma::mat grad_diagPrecMat  = arma::zeros<arma::mat>(400,400);
    grad_diagPrecMat.diag()     = precMat_k.diag() / (this->regularization_parameter_diagonal_PrecMat * this->regularization_parameter_diagonal_PrecMat);

	//derivative of exponential transformation
	grad_diagPrecMat *= precMat_k;

    return(-grad_diagPrecMat);
}



/*
 * Calculate the gradient of the regularizer of the diagonal covariance matrix
 * for a component c
 *
 * cholesky decomposition of diagonal matrix: L is a diagonal matrix
 */
//arma::mat LaplacianApproxRegularizer::gradient_choleskyL_comp_reg(int c){
//
//	arma::mat precMat 	= cinverse_covariance.slice(c);
//	arma::mat gradient_regularizer_diagPrecMat(400,400,arma::fill::zeros);
//	gradient_regularizer_diagPrecMat.diag() = arma::pow(precMat.diag(), -3) / (regularization_diagonal_covMat * regularization_diagonal_covMat); //element wise division of lambda matrix by regularization parameter
//	arma::mat grad_choleskyL(400,400, arma::fill::zeros);
//	arma::mat lower_cholesky = ccholesky_L.slice(c);
//
//	//only diagonal elements count for diagonal matrix
//	for (int k =0; k < 400; k++){
//		for (int l = 0; l <= k ; l++){
//			for (int i = 0; i < 400; i++){
//				grad_choleskyL(k,l)  +=   gradient_regularizer_diagPrecMat(i,k) * lower_cholesky(i,l);
//			}
//			grad_choleskyL(k,l) *= 2;
//		}
//	}
//
//	if(debug  > 1) std::cout << "gradient cholesky L matrix ("<< c << "):\tmin: " << arma::min(arma::vectorise(grad_choleskyL)) <<  ", max: " << arma::max(arma::vectorise(grad_choleskyL)) << ", median:  " << arma::median(arma::vectorise(grad_choleskyL)) << std::endl;
//
//	return(grad_choleskyL);
//}



