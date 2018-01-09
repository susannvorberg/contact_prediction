//============================================================================
// Name        : Parameters.cpp
// Author      : Susann Vorberg
// Version     :
// Copyright   : Your copyright notice
// Description : implementation of parameter class
//============================================================================

#include "Parameters.hpp"

#include <armadillo>
#include <map>

/********************************************************************************************************************************
 *
 * Class constructor
 *
 * @param nr_components_
 *
 *
 *********************************************************************************************************************************/


Parameters::Parameters( int nr_components_)
{

    //set up model specifications
    nr_components	= nr_components_;
    L = 0;
    debug = 0;

    //initialize parameters
    mweight.zeros(2, nr_components);
    mmean.zeros(400, nr_components);
    precMat.zeros(400,400, nr_components);
    covMat.zeros(400,400, nr_components);
    log_det_inverse_covariance.zeros(nr_components);

}

/********************************************************************************************************************************
 *
 * Default Class constructor
 *
 *********************************************************************************************************************************/

Parameters::Parameters()
{
}



/*
 * Print parameters
 */
void Parameters::print(int contact, int component){

	double weight 		    = this->mweight(contact, component);
	arma::vec mean		    = this->mmean.col(component);
	arma::mat covariance	= this->covMat.slice(component);


	std::cout << std::endl << "component " << component << ", contact: " <<  contact << std::endl;
	std::cout << "weight: " << weight << std::endl;
	std::cout << "mean vector, min: "<< arma::min(mean) <<  " ,max: " << arma::max(mean) << ", median:  " << arma::median(mean) << std::endl;
	std::cout << "cov Matrix, min:  "<< arma::min(arma::vectorise(covariance)) <<  ", max: " << arma::max(arma::vectorise(covariance)) << ", median:  " << arma::median(arma::vectorise(arma::trimatu(covariance))) << std::endl << std::endl;

}

/*
 *  Return the number of Gaussian components in the mixture model
 */
int Parameters::get_nr_components(){
    return(this->nr_components);
}

/*
 *  Return the weight of a component either for contacts or non-contacts
 */
double Parameters::get_weight(int component, int contact){
    return(this->mweight(contact, component));
}

/*
 *  Return the mean of a component
 */
arma::vec Parameters::get_mean(int component){
    return(this->mmean.col(component));
}

/*
 *  Return the pecision matrix of a component
 */
arma::mat Parameters::get_precMat(int component){
    return(this->precMat.slice(component));
}

/*
 *  Return the covariance matrix of a component
 */
arma::mat Parameters::get_covMat(int component){
    return(this->covMat.slice(component));
}

/*
 *  Return the logarithm of the determinant of the covariance matrix of a component
 */
double Parameters::get_log_det_inv_covMat(int component){
    return(this->log_det_inverse_covariance(component));
}

/*
 * Define the level of debugging information
 */
void Parameters::set_L(int L){
    this->L = L;
}

/*
 * Define the level of debugging information
 */
void Parameters::set_debug_mode(int debug_mode){
    this->debug = debug_mode;
}

/*
 * Set the parameters that will be optimized
 *
 * mweight 	-> component weights: positive and sum up to 1 for each class (contact, bg)
 * 				ensured by representing them as softmax function
 *
 * mmean	-> mean of gaussian component
 *
 * prec -> precision matrix
 * 			   diagonal == var(X,X) == positive
 * 			   ensured via exp()
 * 			   either isotrope, diagonal or full
 *
 */
void Parameters::set_parameters(std::map<std::string, std::vector<double> > &parameterMap){

	if(debug > 0) std::cout << "Initialise parameters..." << std::endl;


	typedef std::map<std::string, std::vector<double> >::iterator it_type;
	for(it_type iterator = parameterMap.begin(); iterator != parameterMap.end(); iterator++) {

		std::string parameter_name = iterator->first;
		std::vector<double> parameter = iterator->second;
		int parameter_len = parameter.size();

		//determine component number
		int position = parameter_name.rfind ('_');
		int component = std::stoi(parameter_name.substr(position+1));


		if(component >= this->nr_components){
			std::cout << "Parameter name for parameter (" << parameter_name << ") is not valid for model with " << this->nr_components << " components: "  << component << std::endl;
			return;
		}

		if(debug > 0){
		    std::cout << "Parameter: " <<  parameter_name << std::endl;
		    std::cout << "Parameter size: " <<  parameter_len << std::endl;
		    std::cout << "component: " <<  component << std::endl;
		}

		//make sure that parameter values have correct length
		if (parameter_name.find("weight_bg") != std::string::npos) {
			if(parameter_len  != 1){
				std::cout << "Parameter List for parameter " << parameter_name << "is not of length 1 ! It is of length"<< parameter_len  << std::endl;
				return;
			}
			this->mweight(0,component) = exp(parameter[0]);
		}
		else if (parameter_name.find("weight_contact") != std::string::npos) {
			if(parameter_len  != 1){
				std::cout << "Parameter List for parameter " << parameter_name << "is not of length 1 ! It is of length"<< parameter_len << std::endl;
				return;
			}
			this->mweight(1,component) = exp(parameter[0]);
		}
		else if(parameter_name.find("mu") != std::string::npos){
			if(parameter_len  != 400){
				std::cout << "Parameter List for parameter " << parameter_name << "is not of length 400 ! It is of length"<< parameter_len << std::endl;
				return;
			}
			this->mmean.col(component) = arma::vec(parameter);
		}
		else if(parameter_name.find("prec") != std::string::npos){
			if((parameter_len  != 400)){
				std::cout << "Parameter List for parameter " << parameter_name << " is not of length 400 "<< parameter_len << std::endl;
				return;
			}


            //in case of defining lambda_w depending on L
            if (this->L > 0 ){
                //std::cout << "lambda_w determined wrt protein length L= " << this->L << std::endl;
                //std::cout << "parameter " << arma::min(arma::vec(parameter)) << " " << arma::max(arma::vec(parameter)) << std::endl;
                this->precMat.slice(component)  = arma::diagmat( arma::exp(arma::vec(parameter)) * this->L  );
            }else{
                //std::cout << "lambda_w is NOT determined wrt protein length L " << this->L << std::endl;
                this->precMat.slice(component)  = arma::diagmat( arma::exp(arma::vec(parameter)) );
            }

            this->covMat.slice(component) 	= arma::diagmat(1.0 / this->precMat.slice(component).diag());

            //std::cout << "log det " << std::endl;
			//determine the determinant of lambda_k
			this->log_det_inverse_covariance(component) = arma::sum(arma::log(this->precMat.slice(component).diag()));

		}
		else{
			std::cout << "Parameter name does not fit: " << parameter_name << std::endl;
			return;
		}


	}//end for

    //std::cout << "post processing " << std::endl;
	double sum_expweight_bg 		= 0;
	double sum_expweight_contact 	= 0;
	for(int c = 0; c < this->nr_components; c++){
		sum_expweight_bg 		+= this->mweight(0,c); //sum of exp(weights) for background components, weight 1 is always exp(1)!!!
		sum_expweight_contact 	+= this->mweight(1,c); //sum of exp(weights) for contact components, weight 1 is always exp(1)!!!
	}


	//in order to ensure that the weigths are in the range [0-1] and that weights for contact (resp non-contacts) sum to one, we represent them as softmax:
	for(int c = 0; c < this->nr_components; c++){
		this->mweight(0,c) 	/= sum_expweight_bg; 				//exp(weight) / sum of exp(weights) for background components
		this->mweight(1,c) 	/= sum_expweight_contact; 			//exp(weight) / sum of exp(weights) for contact components
	}

	if(debug > 0){

		for(int c = 0; c < this->nr_components; c++){
			this->print(1,c);
			this->print(0,c);
		}
	}

}


