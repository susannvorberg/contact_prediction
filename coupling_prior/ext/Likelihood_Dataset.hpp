//============================================================================
// Name        : Likelihood_Dataset.hpp
// Author      : susi
// Version     :
// Copyright   : Your copyright notice
// Description : return neg ll and gradients for dataset of proteins
//============================================================================

#ifndef LIKELIHOOD_DATASET
#define LIKELIHOOD_DATASET


#include <boost/python.hpp>
#include <armadillo>
#include <map>



//header for class
class Likelihood_Dataset{

	public:
		//constructor
		Likelihood_Dataset(	boost::python::dict pairs_per_protein_,
						    boost::python::dict parameters_
						);


		/*
		 * return the neg LL of the protein dataset
		 * according to specified parameters
		 */
		double get_f();


		/*
         * return the gradients for all parameters in the likelihood function
         * according to order of parameter names list
         */
        boost::python::dict get_gradient_dict();

        /*
         *  Actual computation of neg loglikelihood
         *  NO computation of gradients
         */
        void compute_f();

		/*
		 * Actual computation of neg loglikelihood and gradients for whole dataset
		 */
		void compute_f_df(int hessian_pseudocount);

        /*
        * Set the level of printing debugging information
        */
        void set_debug_mode(int debug_mode);


        /*
        * Set the number of threads for parallelization
        */
        void set_threads_per_protein(int threads);

	private:

		int nr_proteins;
		int nr_components;
		int debug_mode;
		int nr_threads_prot;

		struct myProtein {
		  std::string name;
		  int N;
		  int L;
		  std::string brawfilename;
		  std::string qijabfilename;
		  arma::uvec residue_i;
          arma::uvec residue_j;
          arma::uvec contact;
		} ;

		std::vector<myProtein> dataset;

		double f;
		arma::vec grad_weight_bg;
		arma::vec grad_weight_contact;
		arma::mat grad_mu;
		arma::mat grad_precMat;

		std::map<std::string, std::vector<double> > parameterMap;
		std::map<std::string, std::vector<double> > gradientMap;

		/*
		 * Transform the dictionary/list from python
         * into a map of std::vec of protein structs
		 */
		void protein_dict_to_protein_map(boost::python::dict &pairs_per_protein);

};










#endif //LIKELIHOOD_DATASET

