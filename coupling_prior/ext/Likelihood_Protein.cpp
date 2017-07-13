//============================================================================
// Name        : Likelihood_Protein.cpp
// Author      : Susann Vorberg
// Version     :
// Copyright   : Your copyright notice
// Description :Compute LL and gradients for one protein
//============================================================================


#include "Likelihood_Protein.hpp"
#include "Parameters.hpp"

#define ARMA_64BIT_WORD
#include <armadillo>
#include <map>
#include <msgpack.hpp>

#include <string>
#include <vector>
#include <exception>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>



/********************************************************************************************************************************
 *
 * Class constructor
 *
 *********************************************************************************************************************************/

Likelihood_Protein::Likelihood_Protein(
                    std::string protein_id_,
                    std::string brawfilename_,
                    std::string qijabfilename_,
                    std::map<std::string, std::vector<double> > parameterMap,
                    bool L_dependent
                   )
:protein_id(protein_id_),
brawfilename(brawfilename_),
qijabfilename(qijabfilename_),
parameters(parameterMap.size()/ 4),
L_dependent(L_dependent)
{

    //set up model specifications
    this->debug			    = 0;
    this->nr_components     = parameters.get_nr_components();

    try	{
        read_braw();
        arma::vec diag_regularizer_w 		= arma::vec(400);
        this->regularizer_w =  arma::diagmat(diag_regularizer_w.fill(lambda_w));
    }
    catch (std::exception& e){
        std::cout << "Error reading in braw ("<< protein_id <<") : " << e.what() <<  std::endl;
        exit (EXIT_FAILURE);
    }
    //std::cout << "finished braw..."  <<  std::endl;


    try {set_qijab();}	//sets qijab and Nij
    catch (std::exception& e){
        std::cout << "Error reading in qijab ("<< protein_id <<"): " << e.what() <<  std::endl;
        exit (EXIT_FAILURE);
    }
    //std::cout << "finished qij..."  <<  std::endl;

    //determine precision dependent on protein length L
    if (L_dependent) {
        if (this->debug) std::cout << "Frecision depends on protein length."  <<  std::endl;
        parameters.set_L(this->L);
    }

    //initialize the hyperparameters of the model
    parameters.set_parameters(parameterMap);

    if (this->debug) std::cout << "Finished initialising parameters."  <<  std::endl;

}


/********************************************************************************************************************************
 *
 * Default Class constructor
 *
 *********************************************************************************************************************************/

Likelihood_Protein::Likelihood_Protein()
:protein_id(0),
brawfilename(0),
qijabfilename(0),
L_dependent(0)
{
	std::cout << "Default constructor... " << std::endl;
}


/********************************************************************************************************************************
 *
 * Getter Functions for various controlling features
 *
 *********************************************************************************************************************************/


 /*
 * Get protein length
 */
int Likelihood_Protein::get_L() const{
	return(this->L);
}


/*
 * Return a coupling value for residues i and j and a*20+b
 */
double Likelihood_Protein::get_wijab(int i, int j, int ab) const{

	return(this->w_ij3d(i,j,ab));
}

/*
* Get indices for first residue
*/
arma::uvec Likelihood_Protein::get_i_indices() const{
    return this->i_indices;
}

/*
* Get indices for second residue
*/
arma::uvec Likelihood_Protein::get_j_indices() const{
    return this->j_indices;
}


/*
* Get contact information for pairs
*/
arma::uvec Likelihood_Protein::get_pair_contact_info() const{
    return this->protein_contacts;
}


/*
 * Return the vector q_ij for a given pair i,j
 * see eq 36.: 	qijab := q(x_i = a, x_j = b ) - lambda_w * w_ijab / N_ij
 */
arma::vec Likelihood_Protein::get_qij(int i, int j) const{
	//arma::vec q_ij = q_ij3d.tube(i,j);
	//return q_ij;
	int lin_index = i*(L - (i+1)/2.0 - 1) + j - 1; //i*L - i*(i+1)/2 + j-(i+1);
	arma::vec qij = this->mq_ij.col(lin_index);
	return(qij);
}

/*
 * Return the number_effective sequences for pair i,j
 *
 * --> sum of sequence weights for all sequences that
 * do not contain a gap at position i OR at position j
 * see eq 13: sum_1_N: I(x_i != 0, x_j != 0)
 */
double Likelihood_Protein::get_N_i_j(int i, int j) const {
	return(this->mN_ij(i,j));
}


/*
 * Get values of the neg log likelihood for
 * all pairs in the batch
 */
arma::vec Likelihood_Protein::get_neg_log_likelihood_pairwise() const{

	return(-log_likelihood);
}



/*
 * Get values of the neg log likelihood for the
 * current batch (either only contacting pairs or non-contacting)
 * of the current protein
 */
double Likelihood_Protein::get_neg_log_likelihood() const{

    double neg_log_likelihood_allpairs = -arma::sum(log_likelihood);
	return(neg_log_likelihood_allpairs);
}

/*
 * Return the responsibilities p(k | ij) for all components k for the current pair
 * see eq 67
 * need to run compute_f_df first!
 */
arma::vec Likelihood_Protein::get_responsibilities() const{

	return responsibilities;
}



/*
 * Return gradient of weight
 *  for specified component and contact/noncontact
 *  gradient is cumulated over all pairs
 */
double Likelihood_Protein::get_gradient_weight_comp(int component, int contact){

    //gradient of NEGATIVE loglikelihood
    return -this->grad_weight(component, contact);

}

/*
 * Return gradient of mean vector
 *  for specified component
 *  gradient is cumulated over all pairs
 */
arma::vec Likelihood_Protein::get_gradient_mu_comp(int component){

    //gradient of NEGATIVE loglikelihood
    arma::vec vec_grad_mu = this->grad_mu.col(component);
    return -vec_grad_mu;
}

 /*
 * Return gradient of precision matrix (for NEG log likelihood)
 *  for specified component
 *  gradient is cumulated over all pairs

 * if L_dependent:
 * precision matrix is isotrop and dependending on protein length L
 * e.g. when diag(precMat) == 0.2 * (L-1)
 *
 */
arma::vec Likelihood_Protein::get_gradient_precisionMatrix_comp(int component){

    arma::vec vec_grad_precMat = this->grad_precMat.col(component);

    if (this->L_dependent){
        vec_grad_precMat *= this->L;
    }

    //gradient of NEGATIVE loglikelihood
    return -vec_grad_precMat;
}


/********************************************************************************************************************************
 *
 * Setter functions
 *
 *********************************************************************************************************************************/


/*
* Set the level of printing debugging information
*/
void Likelihood_Protein::set_debug_mode(int debug_mode)
{
	this->debug = debug_mode;
}


/*
 * Set the batch residue pairs
 */
void Likelihood_Protein::set_pairs(arma::uvec &i_indices_in, arma::uvec &j_indices_in, arma::uvec &contact_in)
{

    if (i_indices_in.n_elem == 0 or  j_indices_in.n_elem ==0 or contact_in.n_elem == 0){
        std::cout << "set_pairs: no elements in either i_indices_in, j_indices_in or contact_in for protein " << protein_id << std::endl;
        std::cout << "i_indices_in: " << i_indices_in.n_elem << std::endl;
        std::cout << "j_indices_in: " << j_indices_in.n_elem << std::endl;
        std::cout << "contact_in: " << contact_in.n_elem << std::endl;
    }

    if(i_indices_in.n_elem != j_indices_in.n_elem or i_indices_in.n_elem != contact_in.n_elem){
        throw "Likelihood_Protein::set_pairs: number of indices for i,j or contacts does not match!";
    }

	int max_i = arma::max(i_indices_in);
	int max_j = arma::max(j_indices_in);


	if(max_i >= this->L or max_j >= this->L ){
		throw "Likelihood_Protein::set_pairs: specified amino acid indices do not exist (> protein length L).";
	}

	if (i_indices_in.n_elem != j_indices_in.n_elem or i_indices_in.n_elem != contact_in.n_elem){
	    throw "Likelihood_Protein::set_pairs: Length of indices does not match or is not equal to size of contact vector.";
	}

	this->i_indices = i_indices_in;
	this->j_indices = j_indices_in;
	this->protein_contacts    = contact_in;
	this->number_of_pairs = i_indices.n_elem;
}

/********************************************************************************************************************************
 *
 * Read in couplings and model probabilities
 *
 *********************************************************************************************************************************/

/*
 * Read in braw file and extract all w_ij[1:20, 1:20]
 */
void  Likelihood_Protein::read_braw() {

	if (this->debug) std::cout << "Read in braw File... " << std::endl;

	//deserialize
	msgpack::unpacked msg;

	//find out whether braw file is zipped or not
	if (brawfilename.find("gz") != std::string::npos) {

		if (this->debug) std::cout << "braw  "<< brawfilename <<" compressed!" << std::endl;

		std::ifstream file;
		std::stringstream in;

		//try to open file as fstream
		file.open(brawfilename.c_str(), std::ios_base::in | std::ios_base::binary);

		if (!file) {
			std::cout << "Cannot open gzipped wijab file " + brawfilename  << std::endl;
			exit (EXIT_FAILURE);
		}


		try {
			boost::iostreams::filtering_istreambuf decompressor;
			decompressor.push(boost::iostreams::gzip_decompressor());
			decompressor.push(file);
			boost::iostreams::copy(decompressor, in);
		}
		catch(const boost::iostreams::gzip_error& e) {
	    	std::cout << "Cannot unzip braw file " + brawfilename << ": " << e.what() << std::endl;
	    	exit (EXIT_FAILURE);
		}

        //close ifstream file handle
		file.close();

		// Get the size of the file in bytes
		long fileSize = in.str().size();

		// Convert strstream to const char*
		const std::string& tmp = in.str();
		const char *buf = tmp.c_str();

		msgpack::unpack(&msg, buf, fileSize, 0);

	}else{
		if (this->debug) std::cout << "braw is uncompressed!" << std::endl;
		exit (EXIT_FAILURE);
    }


    //msgpack::object obj = msg.get();
    //std::cout << "obj type: " << obj.type << std::endl;
    //std::cout << "map type: " << msgpack::type::MAP << std::endl;
    //std::cout << "obj size: " << obj.via.map.size << std::endl;

    //convert as a map of msgpack objects
    std::map<std::string, msgpack::object> msgpack_map;
    msg.get().convert(&msgpack_map);


    //protein length
    msgpack_map.at("ncol").convert(&this->L);


    //read lambda_value from: braw.meta['workflow'][0]['regularization']['lambda_pair']
    std::map<std::string, msgpack::object> meta;
    msgpack_map.at("meta").convert(&meta);
    msgpack::object workflow = meta.at("workflow");

    std::vector<msgpack::object> workflow_list;
    workflow.convert(&workflow_list);

    std::map<std::string, msgpack::object> first_element;
    workflow_list[0].convert(&first_element);

    if (first_element.count("regularization")) {

        std::map<std::string, msgpack::object> regularization;
        first_element.at("regularization").convert(&regularization);

        regularization.at("lambda_pair").convert(&this->lambda_w);

    } else {
        std::map<std::string, msgpack::object> parameter;
        first_element.at("parameters").convert(&parameter);

        std::map<std::string, msgpack::object> regularization;
        parameter.at("regularization").convert(&regularization);

        regularization.at("lambda_pair").convert(&this->lambda_w);

    }

    if (this->debug) std::cout << "lambda_pair: " << this->lambda_w << std::endl;

    //get msgpack object for "x_pair"
    // the x_pair msgpack object is of type MAP and has keys "i/j" with j>i
    // this amounts to ncol*(ncol-1)/2 values in this map
    msgpack::object xpair_obj = msgpack_map.at("x_pair");

    //initialise cube for couplings
    w_ij3d.zeros(this->L, this->L, 400);

    //iterate over these i/j pairs and extract the next (and last) map of msgpack objects
    //with keys "i", "j" and "x"
    msgpack::object_kv* p(xpair_obj.via.map.ptr);
    for(msgpack::object_kv* const pend(xpair_obj.via.map.ptr + xpair_obj.via.map.size); p < pend; ++p) {

        //std::string k = p->key.as<std::string>();
        //std::cout << "key: " << k << std::endl;
        std::map<std::string, msgpack::object> xpair_ij_map;
        p->val.convert(&xpair_ij_map);

        int i;
        int j;
        std::vector<double> x_vec;

        xpair_ij_map.at("i").convert(&i);
        xpair_ij_map.at("j").convert(&j);
        xpair_ij_map.at("x").convert(&x_vec); //has been written rowwise

        //read columnwise as matrix
        arma::mat x_mat(arma::vec(x_vec).begin(),21,21);    //length = 21, column-wise to matrix
        arma::mat x_mat_sub = x_mat.submat(0, 0, 19, 19);   //first column 0, last column 19

        //write column-wise back to vector
        w_ij3d.tube(i,j) = arma::vectorise(x_mat_sub);

    }
}

/*
 * Read in q_ijab and Nij from qijab file
 *
 * Necessary to compute Hessian
 */
void Likelihood_Protein::set_qijab() {

	std::ifstream file;
	std::stringstream in;

	//try to open file as ifstream
	file.open(qijabfilename.c_str(), std::ios_base::in | std::ios_base::binary);


	if (!file) {
		std::cout << "Cannot open gzipped qijab file " + qijabfilename  << std::endl;
		exit (EXIT_FAILURE);
	}

    //use boost::decompressor to read from gzip file
	try {
		boost::iostreams::filtering_istreambuf decompressor;
		decompressor.push(boost::iostreams::gzip_decompressor());
		decompressor.push(file);
		boost::iostreams::copy(decompressor, in);
	}
    catch(const boost::iostreams::gzip_error& e) {
    	std::cout << "Cannot unzip qijab file " + qijabfilename << ": " << e.what() << std::endl;
    	exit (EXIT_FAILURE);
    }

    //close ifstream file handle
	file.close();

	// Get the size of the file in bytes
	long fileSize = in.str().size();

	// Convert strstream to const char*
	const std::string& tmp = in.str();
	const char *buf = tmp.c_str();

	//start deserializing
	msgpack::unpacked msg;
	msgpack::unpack(&msg, buf, fileSize, 0);
	std::map<std::string, msgpack::object> msgpack_map;
	msg.get().convert(&msgpack_map);

    //access of vector at map key "N_ij"
    //look script "compute_alignment_features" how vec Nij has been built up
	std::vector<double> Nij_std;
    msgpack_map.at("N_ij").convert(&Nij_std);
    arma::vec Nij_vec = arma::conv_to< arma::vec >::from(Nij_std);

    //has been written as lower triangular matrix row-wise ==
    //fill upper triangle of matrix (without diagonal) column-wise (--> find works columns-wise)
    mN_ij.ones(L, L);
    mN_ij.diag().zeros();
    arma::uvec indices = arma::find(arma::trimatu(mN_ij) > 0);
    mN_ij.elem(indices) = Nij_vec;


    //access of vector at key "q_ij"
    //qij have been written row-wise for ij
    msgpack::object q_ij_array = msgpack_map.at("q_ij");
    std::vector<double> qij_vec_std;
    q_ij_array.convert(&qij_vec_std);

    mq_ij = arma::conv_to<arma::mat>::from(qij_vec_std);
    mq_ij.reshape(400, L*(L-1)/2);




}





/********************************************************************************************************************************
 *
 * Compute quantities
 *
 *********************************************************************************************************************************/

/*
 * Compute the Hessian for a certain pair ij
 */
arma::mat Likelihood_Protein::compute_H_ij(int i, int j, int contact )
{


    int lin_index       = i*(L - (i+1)/2.0 - 1) + j - 1; //i*L - i*(i+1)/2 + j-(i+1);


    arma::vec vqij = mq_ij.col(lin_index);
    double N_ij = mN_ij(i,j);
    arma::vec w_ij = arma::vec(w_ij3d.tube(i,j));

    arma::mat Qij 		= arma::diagmat(vqij);

    //outer product
    arma::mat qij_prod = vqij * vqij.t();

    //determine negative Hessian = Hij
    arma::mat diff = Qij - qij_prod;
    arma::mat H_ij = N_ij * diff + regularizer_w; //eq 37

    return H_ij;
}

 /*
 * Compute the ratio of the two Gaussians in log space
 *       N(0 | mu_k, lambda_k)
 * ------------------------------
 *   N(0 | mu_ij_k, lambda_ij,k)
 * = - 0.5 [ (log|Sigma1|  + mu1^TSigma1^-1mu1) - ( log|Sigma2|  + mu2^TSigma2^-1mu2 ) ]
 * =   0.5 [ (log|Lambda1| - mu1^T Lambda1 mu1) - ( log|Lambda2| - mu2^T Lambda2 mu2 ) ]
 *
 */
double Likelihood_Protein::log_density_gaussian_ratio(	const arma::vec &mu_k,
                                                        const arma::vec &mu_ij_k,
                                                        const arma::mat &lambda_k,
                                                        const arma::mat &lambda_ij_k,
                                                        double log_det_lambdak,
                                                        double log_det_lambdaijk)
{

	double gaussian_1 = log_det_lambdak - arma::as_scalar(mu_k.t()*lambda_k*mu_k);
	double gaussian_2 = log_det_lambdaijk - arma::as_scalar(mu_ij_k.t()*lambda_ij_k*mu_ij_k);
	double gaussian_ratio_log_density = 0.5 * (gaussian_1 - gaussian_2);


//    if(gaussian_ratio_log_density>1000){
//        std::cout  << "gaussian_1: " << 0.5 * gaussian_1 << std::endl;
//        std::cout  << "gaussian_2: " << 0.5 * gaussian_2 << std::endl;
//    }


/*	if(! std::isnormal(gaussian_1)) {
	    std::cout  << "gaussian_1: " << gaussian_1 << std::endl;
	    std::cout  << "log_det_lambdak: " << log_det_lambdak << std::endl;
	    std::cout  << "arma::as_scalar(mu_k.t()*lambda_k*mu_k): " << arma::as_scalar(mu_k.t()*lambda_k*mu_k) << std::endl;
	}

	if(! std::isnormal(gaussian_2)) {
	    std::cout  << "gaussian_2: " << gaussian_2 << std::endl;
	    std::cout  << "log_det_lambdaijk: " << log_det_lambdaijk << std::endl;
	    std::cout  << "arma::as_scalar(mu_ij_k.t()*lambda_ij_k*mu_ij_k): " << arma::as_scalar(mu_ij_k.t()*lambda_ij_k*mu_ij_k) << std::endl;
	}*/


	return(gaussian_ratio_log_density);
}



/*
 * This function calculates the Hessian and mu_ij_k and Lambda_ij_k
 * for all pairs (i,j) in the batch
 * and for all components k
 *
*/
void Likelihood_Protein::compute_negLL(int nr_threads_prot)
{
    //initialize likelihood vector
    log_likelihood.zeros(this->number_of_pairs);

	#pragma omp parallel for num_threads(nr_threads_prot)
    for(int pair = 0; pair < this->number_of_pairs; pair++){

		int i 				= this->i_indices(pair);
		int j 				= this->j_indices(pair);
		int lin_index       = i*(L - (i+1)/2.0 - 1) + j - 1; //i*L - i*(i+1)/2 + j-(i+1);
		int contact         = this->protein_contacts(pair);

        //for computation of likelihood
		arma::vec log_density(this->nr_components, arma::fill::zeros);

        arma::vec vqij = mq_ij.col(lin_index);
        double N_ij = mN_ij(i,j);
		arma::vec w_ij = w_ij3d.tube(i,j);
		//arma::vec vqij = q_ij3d.tube(i,j);

		//diagonal matrix Qij = diag(q'ji)
		//q'ijab = q(x_i=a, x_j=b) - (lambda_w * wijab / N_ij) --> has been precomputed
		arma::mat Qij 		= arma::diagmat(vqij);

		//outer product
		arma::mat qij_prod = vqij * vqij.t();

		//determine negative Hessian = Hij
		arma::mat diff = Qij - qij_prod;
		arma::mat H_ij = N_ij * diff + regularizer_w; //eq 37

		//precompute product H_ij * wij
		arma::vec Hij_wij_prod = H_ij * w_ij;


		for(int k = 0; k < this->nr_components; k++){


            //gaussian parameters of coupling prior
            double weight_k 		        = this->parameters.get_weight(k, contact);
			arma::vec mu_k 				    = this->parameters.get_mean(k);
			arma::mat lambda_k 			    = this->parameters.get_precMat(k);


            //---------------- simplify computation in case that lambda_k is diagonal matrix
            // A, A_inv, lambda_k, Qij are diagonal matrices
            arma::vec A = N_ij * vqij + lambda_k.diag();     //represents diagonal matrix
            arma::vec A_inv = 1.0 / A;                    //represents diagonal matrix
            double log_det_A	= arma::sum(arma::log(A));
            arma::vec Ainv_qij_product  = A_inv % vqij;   ////column vector
            double triple_product 	    = arma::sum(vqij % Ainv_qij_product);
            //---------------- matrix computations in case lambda_k is NOT diagonal matrix
			//A is a diagonal matrix, as Qij and lambda_k are diagonal matrices
			//arma::mat A 		= N_ij * Qij + lambda_k;                    //diagonal
			//arma::mat A_inv 	= arma::diagmat(1.0 / A.diag());            //diagonal
			//double log_det_A	= arma::sum(arma::log(A.diag()));
			//arma::vec Ainv_qij_product  = arma::vec(A_inv * vqij);          //400x1 dim matrix
			//double triple_product 	    = arma::as_scalar(vqij.t() * Ainv_qij_product);


			//compute lambda_ij_k_mat
			arma::mat lambda_ij_k_mat  	= H_ij - regularizer_w + lambda_k;


            //debugging: we assume diagonal Hessian ================================================================
//            arma::mat lambda_ij_k_mat_inv(400,400,arma::fill::zeros);
//            lambda_ij_k_mat_inv.diag() = 1.0 / lambda_ij_k_mat.diag();
            //debugging=======================================================================================
			//compute inverse of lambda_ij_k_mat
			//---------------- simplify computation in case that lambda_k is diagonal matrix
			arma::mat lambda_ij_k_mat_inv = arma::diagmat(A_inv) + (Ainv_qij_product * Ainv_qij_product.t()) / (1.0/N_ij - triple_product);
			//---------------- matrix computations in case lambda_k is NOT diagonal matrix
			//arma::mat  lambda_ij_k_mat_inv  = A_inv + (Ainv_qij_product * Ainv_qij_product.t()) / (1.0/N_ij - triple_product);



		    //compute mu_ij_k from precomputed entities
			arma::vec mu_ij_k_vec      = lambda_ij_k_mat_inv * ( Hij_wij_prod + lambda_k * mu_k);



            //debugging: we assume diagonal Hessian ================================================================
//            log_det_lambda_ij_k(k) = arma::sum(arma::log(lambda_ij_k_mat.diag()));
            //debugging=======================================================================================
			//save log determinant of lambda_ij_k, see page 16 Saikats theory
			double log_det_lambda_ij = log(1 - N_ij * triple_product) + log_det_A;

			//ratio of two gaussians in log space
			//     N(0 | mu_k, lambda_k)
            //------------------------------
            //  N(0 | mu_ij_k, lambda_ij,k)
			double gaussian_ratio_logdensity = log_density_gaussian_ratio(	mu_k,
																			mu_ij_k_vec,
																			lambda_k,
																			lambda_ij_k_mat,
																			this->parameters.get_log_det_inv_covMat(k),
																			log_det_lambda_ij);


            log_density(k) = log(weight_k) + gaussian_ratio_logdensity;

		}//end loop over components k



		//Johannes suggestion how to precisely compute the responsibilities
		double a_max = arma::max(log_density);
		arma::vec resps = arma::exp(log_density - a_max);//r_nk = exp( a_nk - a_max)
		double sum_resp = arma::sum(resps);    //sum += r_nk


		//save neg likelihood of current pair
		double f = log(sum_resp) + a_max;

		if(! std::isnormal(f)) {
				std::cout  << "ERROR: likelihood cannot be computed for protein " << protein_id << ", i " << i << ", j " << j << " ("<< contact <<"): " << f << std::endl;
                std::cout  << "Nij " << N_ij << ", sum_resp: " << sum_resp << ", a_max: " << a_max << std::endl;
                for(int k = 0; k < this->nr_components; k++){
                    std::cout  << "component: " << k << ", sum_precMat(k)diag: "<< arma::sum(this->parameters.get_precMat(k).diag()) << ", responsibilty:" << resps(k)/sum_resp << ", log_density: " << log_density(k) << std::endl;
                }

				continue;
		} else log_likelihood(pair) = f;

	}//end of parallelized for loop over ij pairs

}


/*
 * Compute the following entities:
 *      - responsibilities r(k|ij)
 *      - mu_ij_k
 *      - Lambda_ij_k
 *      - log determinante Lamda_ij_k
 *
 *  And then use these to compute likelihood and gradients
 */
void Likelihood_Protein::compute_f_df(int hessian_pseudocount)
{


	//initialize parameters for second gaussian
    mu_ij_k.zeros(400, this->nr_components);
    lambda_ij_k.zeros(400, 400, this->nr_components);
    lambda_ij_k_inv.zeros(400, 400, this->nr_components);
    log_det_lambda_ij_k.zeros(this->nr_components);
	responsibilities.zeros(this->nr_components);

    //initialize likelihood vector
    log_likelihood.zeros(this->number_of_pairs);

    //intialize gradients to zero
	grad_weight.zeros(this->nr_components, 2);
	grad_mu.zeros(400,this-> nr_components);
	grad_precMat.zeros(400, this->nr_components);

    for(int pair = 0; pair < this->number_of_pairs; pair++){

		int i 				= this->i_indices(pair);
		int j 				= this->j_indices(pair);
		int lin_index       = i*(L - (i+1)/2.0 - 1) + j - 1; //i*L - i*(i+1)/2 + j-(i+1);
		int contact         = this->protein_contacts(pair);

        //for computation of likelihood
		arma::vec log_density(this->nr_components, arma::fill::zeros);

        arma::vec vqij = mq_ij.col(lin_index);
        double N_ij = mN_ij(i,j);
		arma::vec w_ij = w_ij3d.tube(i,j);
		//arma::vec vqij = q_ij3d.tube(i,j);

		//diagonal matrix Qij = diag(q'ji)
		//q'ijab = q(x_i=a, x_j=b) - (lambda_w * wijab / N_ij) --> has been precomputed
		arma::mat Qij 		= arma::diagmat(vqij);

		//outer product
		arma::mat qij_prod = vqij * vqij.t();

		//determine negative Hessian = Hij
		arma::mat diff = Qij - qij_prod;
		arma::mat H_ij = N_ij * diff + regularizer_w; //eq 37

        // remove after debugging ---------------------------------------------------------------------------------
		//debugging artefacts in learning
		//add counts to Hessian to see if we have problems with non-informative data
		H_ij.diag() += hessian_pseudocount;
		//H_ij = arma::diagmat(H_ij);
		// remove after debugging ---------------------------------------------------------------------------------

		//precompute product H_ij * wij
		arma::vec Hij_wij_prod = H_ij * w_ij;

		for(int k = 0; k < this->nr_components; k++){

            //gaussian parameters of coupling prior
            double weight_k 		        = this->parameters.get_weight(k, contact);
			arma::vec mu_k 				    = this->parameters.get_mean(k);
			arma::mat lambda_k 			    = this->parameters.get_precMat(k);


            //---------------- simplify computation in case that lambda_k is diagonal matrix
            // A, A_inv, lambda_k, Qij are diagonal matrices
            arma::vec A = N_ij * vqij + lambda_k.diag();     //represents diagonal matrix
            arma::vec A_inv = 1.0 / A;                    //represents diagonal matrix
            double log_det_A	= arma::sum(arma::log(A));
            arma::vec Ainv_qij_product  = A_inv % vqij;   ////column vector
            double triple_product 	    = arma::sum(vqij % Ainv_qij_product);
            //---------------- matrix computations in case lambda_k is NOT diagonal matrix
			//A is a diagonal matrix, as Qij and lambda_k are diagonal matrices
			//arma::mat A 		= N_ij * Qij + lambda_k;                    //diagonal
			//arma::mat A_inv 	= arma::diagmat(1.0 / A.diag());            //diagonal
			//double log_det_A	= arma::sum(arma::log(A.diag()));
			//arma::vec Ainv_qij_product  = arma::vec(A_inv * vqij);          //400x1 dim matrix
			//double triple_product 	    = arma::as_scalar(vqij.t() * Ainv_qij_product);


			//compute lambda_ij_k_mat
			arma::mat lambda_ij_k_mat  	= H_ij - regularizer_w + lambda_k;
			lambda_ij_k.slice(k) 	    = lambda_ij_k_mat;


            //debugging: we assume diagonal Hessian ================================================================
//            arma::mat lambda_ij_k_mat_inv(400,400,arma::fill::zeros);
//            lambda_ij_k_mat_inv.diag() = 1.0 / lambda_ij_k_mat.diag();
            //debugging=======================================================================================
			//compute inverse of lambda_ij_k_mat
			//---------------- simplify computation in case that lambda_k is diagonal matrix
			arma::mat lambda_ij_k_mat_inv = arma::diagmat(A_inv) + (Ainv_qij_product * Ainv_qij_product.t()) / (1.0/N_ij - triple_product);
			//---------------- matrix computations in case lambda_k is NOT diagonal matrix
			//arma::mat  lambda_ij_k_mat_inv  = A_inv + (Ainv_qij_product * Ainv_qij_product.t()) / (1.0/N_ij - triple_product);
			lambda_ij_k_inv.slice(k) 	    = lambda_ij_k_mat_inv;


		    //compute mu_ij_k from precomputed entities
			arma::vec mu_ij_k_vec      = lambda_ij_k_mat_inv * ( Hij_wij_prod + lambda_k * mu_k);
			mu_ij_k.col(k)             = mu_ij_k_vec;


            //debugging: we assume diagonal Hessian ================================================================
//            log_det_lambda_ij_k(k) = arma::sum(arma::log(lambda_ij_k_mat.diag()));
            //debugging=======================================================================================
			//save log determinant of lambda_ij_k, see page 16 Saikats theory
			log_det_lambda_ij_k(k) = log(1 - N_ij * triple_product) + log_det_A;

			//ratio of two gaussians in log space
			//     N(0 | mu_k, lambda_k)
            //------------------------------
            //  N(0 | mu_ij_k, lambda_ij,k)
			double gaussian_ratio_logdensity = log_density_gaussian_ratio(	mu_k,
																			mu_ij_k_vec,
																			lambda_k,
																			lambda_ij_k_mat,
																			this->parameters.get_log_det_inv_covMat(k),
																			log_det_lambda_ij_k(k) );



            if(this->debug > 0 and gaussian_ratio_logdensity>1000){
                std::cout  << protein_id << " " << i << " " << j << " " << pair << " " << contact << " " << k << std::endl;
                std::cout  << "Gaussian log density > 1: " << gaussian_ratio_logdensity << std::endl;
                std::cout  << "A_inv.max() : " << A_inv.max()  << std::endl;
                arma::uword max_ind = index_max(A_inv);
                std::cout  << "A(max_ind) : " << A(max_ind)  << std::endl;
                std::cout  << "vqij(max_ind): " << vqij(max_ind) << std::endl;
                std::cout  << "N_ij: " << N_ij << std::endl;
                std::cout  << "mu_ij_k_vec(max_ind) : " << mu_ij_k_vec(max_ind)  << std::endl;
                std::cout  << "mu_k(max_ind) : " << mu_k(max_ind)  << std::endl;
                std::cout  << "lambda_ij_k_mat(max_ind, max_ind) : " << lambda_ij_k_mat(max_ind, max_ind)  << std::endl;
                std::cout  << "lambda_ij_k_mat_inv(max_ind, max_ind) : " << lambda_ij_k_mat_inv(max_ind, max_ind)  << std::endl;
                std::cout  << "lambda_k(max_ind, max_ind) : " << lambda_k(max_ind, max_ind)  << std::endl;
                std::cout  << "parameters.get_log_det_inv_covMat(k) : " << this->parameters.get_log_det_inv_covMat(k)  << std::endl;
                std::cout  << "log_det_lambda_ij_k(k) : " << log_det_lambda_ij_k(k)  << std::endl;
            }


            log_density(k) = log(weight_k) + gaussian_ratio_logdensity;

		}//end loop over components k

		//Johannes suggestion how to precisely compute the responsibilities
		double a_max = arma::max(log_density);
		this->responsibilities = arma::exp(log_density - a_max);//r_nk = exp( a_nk - a_max)
		double sum_resp = arma::sum(this->responsibilities);    //sum += r_nk
		this->responsibilities /= sum_resp; //normalization of responsibilities, NOT in log space  => r_nk /= sum;



		//save neg likelihood of current pair
		double f = log(sum_resp) + a_max;
		if(! std::isnormal(f)) {
				std::cout  << "ERROR: likelihood cannot be computed for protein " << protein_id << ", i " << i << ", j " << j << " ("<< contact <<"): " << f << std::endl;
                std::cout  << "Nij: " << N_ij << ", sum_resp: " << sum_resp << ", a_max: " << a_max << std::endl;
                for(int k = 0; k < this->nr_components; k++){
                    std::cout  << "component: " << k << ", sum_precMat(k)diag: "<< arma::sum(this->parameters.get_precMat(k).diag()) << ", responsibilty:" << this->responsibilities(k) << ", log_density: " << log_density(k) << ", log_det_lambda_ij_k: " << log_det_lambda_ij_k(k) << std::endl;
                }

				continue;
		} else log_likelihood(pair) = f;




        // Compute ALL gradients for ALL components
        for(int k = 0; k < this->nr_components; k++){

            //weights (either for contact or non_contact)
            grad_weight(k, contact) += gradient_weight_comp(k, contact);

            //mu
            grad_mu.col(k) += gradient_mu_comp(k);

            //precMat - diagonal
            arma::mat grad_precMat_protein = gradient_precisionMatrix_comp(k);
            grad_precMat.col(k) += grad_precMat_protein.diag();

        }//end components


	}//end loop over ij pairs


}






/********************************************************************************************************************************
 *
 * Compute gradients
 *
 *********************************************************************************************************************************/



/*
 * Calculate the gradient of mu
 * for component k and the current residue pair
 *
 * @return: mu vector: 400d
 */
arma::vec Likelihood_Protein::gradient_mu_comp(int k)
{

	//parameters of coupling prior
	arma::mat lambda_k  = this->parameters.get_precMat(k);
	arma::vec mu_k 		= this->parameters.get_mean(k);

    //parameters of second gaussian
    arma::vec mu_ij_k_vec     = this->mu_ij_k.col(k);

    arma::vec diff            = mu_ij_k_vec - mu_k;

    //actual gradient calculation
    //arma::vec grad_mu_k = responsibilities(k) * lambda_k * (mu_ij_k_vec - mu_k);
    //in case of diagonal lambda_k : element-wise multiplication
    arma::vec grad_mu_k = this->responsibilities(k) * (lambda_k.diag() % diff);

    //this is the gradient for the LOG LIKELIHOOD
	return(grad_mu_k);

}


/*
 * This function calculates the gradient for
 * the precision matrix Lambda
 * for component k and current residue pair
 *
 * @return: precision matrix Lambda: 400 x 400
 */
arma::mat Likelihood_Protein::gradient_precisionMatrix_comp(int k)
{

	//parameters of coupling prior
	arma::mat covMat    = this->parameters.get_covMat(k);
	arma::mat lambda_k 	= this->parameters.get_precMat(k);
	arma::vec mu_k 	    = this->parameters.get_mean(k);

    //parameters of second gaussian
    arma::mat lambda_ij_k_mat_inv   = this->lambda_ij_k_inv.slice(k);
    arma::vec mu_ij_k_vec           = this->mu_ij_k.col(k);

    arma::vec diff                  = mu_ij_k_vec - mu_k;

    //outer product of vector with its tranpose
    //arma::mat diffmat               = diff * diff.t();
    //arma::mat grad_precisionMatrix_k  = 0.5 * responsibilities(k) *  (covMat  - lambda_ij_k_mat_inv - diffmat );
    //in case of diagonal precision matrix, only diagonals matter
    arma::mat grad_precisionMatrix_k(400,400,arma::fill::zeros);
    grad_precisionMatrix_k.diag() = 0.5 * this->responsibilities(k) * (covMat.diag() - lambda_ij_k_mat_inv.diag() - arma::square(diff) );

	//derivative of exponential transformation
	//grad_precisionMatrix_k *= lambda_k;
	//in case of diagonal matrix
	grad_precisionMatrix_k.diag() %= lambda_k.diag();



//    double max_val = arma::abs(grad_precisionMatrix_k).max();
//	if(max_val > 1.0){
//        std::cout  << "max gradient for prec Mat > 1 for component " << k << "!" << std::endl;
//        arma::uword i = arma::abs(grad_precisionMatrix_k).index_max();
//        std::cout  << "index of max:  " << i << std::endl;
//        std::cout  << "grad_precisionMatrix_k[a,b]:  " << grad_precisionMatrix_k(i) << std::endl;
//        std::cout  << "max_val:  " << max_val << std::endl;
//        std::cout  << "precMat[a,b]:  " << precMat(i) << std::endl;
//        std::cout  << "responsibilities(k):  " << responsibilities(k) << std::endl;
//        std::cout  << "covMat[a,b]:  " << covMat(i) << std::endl;
//        std::cout  << "lambda_ij_k_mat_inv[a,b]:  " << lambda_ij_k_mat_inv(i) << std::endl;
//        std::cout  << "diffmat[a,b]:  " << diffmat(i) << std::endl;
//	}

    //this is the gradient for LOG LIKELIHOOD
	return(grad_precisionMatrix_k);
}





/*
 * This function calculates the gradient for the component weight
 * for component k and contact/noncontact, e.g. and the current pair
 */
double Likelihood_Protein::gradient_weight_comp(int k, int contact)
{

    //weight_k is softmax of actual weight parameter --> derive softmax function:
	double weight_k 	    = this->parameters.get_weight(k, contact);
    double grad_weight_k    = this->responsibilities(k)- weight_k;

    //returns the gradient for the LOG LIKELIHOOD
	return(grad_weight_k);

}






