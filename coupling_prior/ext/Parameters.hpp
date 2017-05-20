//============================================================================
// Name        : Parameters.hpp
// Author      : susi
// Version     :
// Copyright   : Your copyright notice
// Description : header file for parameters of ll
//============================================================================

#ifndef PARAMETERS
#define PARAMETERS

#include <armadillo>
#include <map>

/********************************************************************************************************************************
 *
 * Class header
 *
 *********************************************************************************************************************************/

class Parameters{
	public:
		//constructor
		Parameters(int nr_components_);

		//default constructor
		Parameters();


        /*
         * Print parameters
         */
        void print(int contact, int component);

        /*
         * Define the level of debugging information
         */
        void set_L(int L);

        /*
         * Define the level of debugging information
         */
        void set_debug_mode(int debug_mode);

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
        void set_parameters(std::map<std::string, std::vector<double> > &parameterMap);

        //getter Functions
        int get_nr_components();
        double get_weight(int component, int contact);
        arma::vec get_mean(int component);
        arma::mat get_precMat(int component);
        arma::mat get_covMat(int component);
        double get_log_det_inv_covMat(int component);

	private:
		int debug;
		int nr_components;
	    int L;
		arma::mat mweight;			//parameters to optimize
		arma::mat mmean;			//parameters to optimize
		arma::cube covMat;		//parameters to optimize
		arma::cube precMat;
		arma::vec log_det_inverse_covariance;
};

#endif //PARAMETERS