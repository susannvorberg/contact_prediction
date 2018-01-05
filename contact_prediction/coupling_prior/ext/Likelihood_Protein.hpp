//============================================================================
// Name        : Likelihood_Protein.hpp
// Author      : Susann Vorberg
// Version     :
// Copyright   : Your copyright notice
// Description : Compute LL and gradients for one protein
//============================================================================


#ifndef LIK_PROTEIN
#define LIK_PROTEIN

#define ARMA_64BIT_WORD


#include "Parameters.hpp"
#include <armadillo>
#include <map>


class Likelihood_Protein{
	public:

	    Likelihood_Protein(	std::string protein_id_,
                            std::string brawfilename_,
                            std::string qijabfilename_,
                            std::map<std::string, std::vector<double> > parameterMap,
                            bool L_dependent
                            );

        //default constructor
		Likelihood_Protein();




    /*
     * Get protein length
     */
    int get_L() const;

    /*
    * Get indices for first residue
    */
    arma::uvec get_i_indices() const;

    /*
    * Get indices for second residue
    */
    arma::uvec get_j_indices() const;

    /*
    * Get contact information for pairs
    */
    arma::uvec get_pair_contact_info() const;

    /*
     * Return a coupling value for residues i and j and a*20+b
     */
    double get_wijab(int i, int j, int ab) const;


    /*
     * Return the vector q_ij for a given pair i,j
     * see eq 36.: 	qijab := q(x_i = a, x_j = b ) - lambda_w * w_ijab / N_ij
     */
    arma::vec get_qij(int i, int j) const;

    /*
     * Return the number_effective sequences for pair i,j
     *
     * --> sum of sequence weights for all sequences that
     * do not contain a gap at position i OR at position j
     * see eq 13: sum_1_N: I(x_i != 0, x_j != 0)
     */
    double get_N_i_j(int i, int j) const;


    /*
     * Get values of the neg log likelihood for
     * all pairs in the batch
     */
    arma::vec get_neg_log_likelihood_pairwise() const;

    /*
     * Get values of the neg log likelihood for the
     * current batch (either only contacting pairs or non-contacting)
     * of the current protein
     */
    double get_neg_log_likelihood() const;

    /*
     * Return the responsibilities p(k | ij) for all components k for the current pair
     * see eq 67
     * need to run compute_f_df first!
     */
    arma::vec get_responsibilities() const;


    /*
     * Return gradient of weight
     *  for specified component and contact/noncontact
     *  gradient is cumulated over all pairs
     */
    double get_gradient_weight_comp(int component, int contact);

    /*
     * Return gradient of mean vector
     *  for specified component
     *  gradient is cumulated over all pairs
     */
    arma::vec get_gradient_mu_comp(int component);

     /*
     * Return gradient of precision matrix (for NEG log likelihood)
     *  for specified component
     *  gradient is cumulated over all pairs

     * if L_dependent:
     * precision matrix is isotrop and dependending on protein length L
     * e.g. when diag(precMat) == 0.2 * (L-1)
     *
     */
    arma::vec get_gradient_precisionMatrix_comp(int component);


    /*
    * Set the level of printing debugging information
    */
    void set_debug_mode(int debug_mode);

    /*
     * Set the indices of residue pairs and their contact information
     */
    void set_pairs(arma::uvec &i_indices_in, arma::uvec &j_indices_in, arma::uvec &contact_in);

    /*
     * Read in braw file and extract all w_ij[1:20, 1:20]
     */
    void read_braw();

    /*
     * Read in q_ijab and Nij from qijab file
     *
     * Necessary to compute Hessian
     */
    void set_qijab();


    /*
     * Compute the Hessian for a certain pair ij
     */
    arma::mat compute_H_ij(int i, int j, int contact );


     /*
     * Compute the ratio of the two Gaussians in log space
     *       N(0 | mu_k, lambda_k)
     * ------------------------------
     *   N(0 | mu_ij_k, lambda_ij,k)
     * = - 0.5 [ (log|Sigma1|  + mu1^TSigma1^-1mu1) - ( log|Sigma2|  + mu2^TSigma2^-1mu2 ) ]
     * =   0.5 [ (log|Lambda1| - mu1^T Lambda1 mu1) - ( log|Lambda2| - mu2^T Lambda2 mu2 ) ]
     *
     */
    double log_density_gaussian_ratio(	const arma::vec &mu_k,
                                        const arma::vec &mu_ij_k,
                                        const arma::mat &lambda_k,
                                        const arma::mat &lambda_ij_k,
                                        double log_det_lambdak,
                                        double log_det_lambdaijk);

    /*
     * This function calculates the Hessian and mu_ij_k and Lambda_ij_k
     * for all pairs (i,j) in the batch
     * and for all components k
     *
    */
    void compute_negLL(int nr_threads_prot);


    /*
     * Compute the following entities:
     *      - responsibilities r(k|ij)
     *      - mu_ij_k
     *      - Lambda_ij_k
     *      - log determinante Lamda_ij_k
     */
    void compute_f_df(int hessian_pseudocount);


    /*
     * Calculate the gradient of mu
     * for component k and the current residue pair
     */
    arma::vec gradient_mu_comp(int k);

    /*
     * This function calculates the gradient for
     * the precision matrix Lambda
     * for component k and current residue pair
     *
     * @return: precision matrix Lambda: 400 x 400
     */
    arma::mat gradient_precisionMatrix_comp(int k);

    /*
     * This function calculates the gradient for the component weight
     * for component k and contact/noncontact, e.g. and the current pair
     */
    double gradient_weight_comp(int k, int contact);



    //protected members can be accesssed by derived class
    private:
		const std::string protein_id;
		const std::string brawfilename;
		const std::string qijabfilename;
    	Parameters parameters;
    	const bool L_dependent;

        arma::vec responsibilities;
        arma::vec log_likelihood;

        arma::cube lambda_ij_k;
        arma::cube lambda_ij_k_inv;
        arma::mat mu_ij_k;
        arma::vec log_det_lambda_ij_k;

        arma::mat grad_weight;
        arma::mat grad_mu;
        arma::mat grad_precMat;

        arma::uvec i_indices;
		arma::uvec j_indices;
		arma::uvec protein_contacts;
		arma::mat mq_ij;
		arma::mat mN_ij;
		arma::cube w_ij3d;
		arma::mat regularizer_w;

        int L;
    	double lambda_w;
    	int number_of_pairs;
		int debug;
		int nr_components;


};



#endif // LIK_PROTEIN



