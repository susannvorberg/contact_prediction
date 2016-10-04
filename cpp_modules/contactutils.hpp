//============================================================================
// Name        : contactutils.hpp
// Author      : susi
// Copyright   : Your copyright notice
// Description : Header file for
//				 Functions needed in python
//				 Made accessable through boost::python
//============================================================================

#ifndef CONTACTUTILS
#define CONTACTUTILS

#include <boost/python.hpp>
#include <armadillo>




/********************************************************************************************************************************
 *
 * General Fct
 *
 *********************************************************************************************************************************/

#define SEQ_BUFFER 8192
const double PI = std::atan(1.0)*4;
const double log2pi = std::log(2.0 * PI);


const unsigned char AMINO_INDICES[26] = {
		     //  A   B   C   D   E   F   G   H   I   J   K   L   M   N   O   P   Q   R   S   T   U   V   W   X   Y   Z
			     1,  0,  5,  4,  7, 14,  8,  9,  10, 0, 12, 11, 13,  3,  0, 15,  6,  2, 16, 17,  0, 20, 18,  0, 19,  0 };

const unsigned char CHAR_INDICES[21] = {
		     //  0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20
		     //  -   A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V
				45, 65, 82, 78, 68, 67, 81, 69, 71, 72, 73, 76, 75, 77, 70, 80, 83, 84, 87, 89, 86};



/*
 * Convert amino acid ASCII code to an index, with 0 = gap and 1 = A, 2 = C, etc.
 */
unsigned char aatoi(unsigned char aa);
unsigned char itoaa(unsigned char i);


/*
 * Given index i and index j of an nxn matrix
 * return the index ij of the vectorized upper triangular matrix
 * without diagonal !!!
 * row-wise !!!
 */
double get_linear_index(int i, int j, int n);

/*
 * From a linear index, rowwise, no diagonal
 * reconstruct paired index
 * nxn Matrix
 */
arma::vec get_paired_index(int ij, int n);

/*
 * Get the size of a file
 */
long getFileSize(FILE *file);



/*
 * Setup  a matrix of size n x n
 * and fill the upper triangle with values from input vector
 * values needs to be of size n*(n+1)/2
 */
arma::mat fill_matrix_from_triangular_indices(arma::vec values, int n, bool indices_upper);





/********************************************************************************************************************************
 *
 * Gaussian Mixtures
 *
 *********************************************************************************************************************************/
/*
 * inefficient way of calculating gaussian density with covMatrix as input
*/
double mydmvnorm_cov_logspace(const arma::vec &x, const arma::vec &meanVector, const arma::mat &cov, bool logd );

/*
 * Fast Implementation of multivariate gaussian density after baysm R package
 * using Prec Matrix as input
 */
arma::vec mydmvnorm_prec_logspace_chol(const arma::mat &x, const arma::rowvec &mean, const arma::mat &prec, bool logd, int cores );

/*
 * Fast Implementation of multivariate gaussian density after baysm R package
 * using Cov Matrix as input
 */
arma::vec dmvnrm_arma_mc(const arma::mat &x, const arma::rowvec &mean, const arma::mat &sigma, bool logd, int cores);



/********************************************************************************************************************************
 *
 * Bayes Functions
 *
 *********************************************************************************************************************************/



/** Read an MSA file into an index matrix.
 *
 * Matrix will be returned as a pointer to a one-dimensional int array
 * where the cell at row i and column j will have the index i + nrow*j
 */
arma::mat read_msa(std::string msafilename, int N, int L , int debug);

/*
 * Read vector of sequence weights from file
 *
 * each element contains the weight of a sequence from the alignment
 */
arma::vec read_sequenceweights(std::string seqweightsfilename, int N, int debug);


/*
 * Read in qij from qijab file for specified i and j
 *
 * qijab and Nij are necessary to compute Hessian
 */
std::vector<double> read_qij_msgpack(std::string qijabfilename, int i, int j, int debug);


/*
 * Read in q_ijab and Nij from qijab file
 *
 * qijab and Nij are necessary to compute Hessian
 */
arma::mat read_Nij_msgpack(std::string qijabfilename, int L, int debug);


/*
 * Read in braw file and extract all w_ij[1:20, 1:20]
 */
arma::cube  read_braw(std::string brawfilename, int L, int debug);



/*
 * Calculate the weights for each sequence based on pairwise amino acid comparisons
 */
arma::vec calculate_seq_weights(int N,int L, arma::mat psicov, double reweighting_threshold, int debug);


/*
 * Calculate the APC correction for a given matrix
 */
arma::mat calc_apc(arma::mat Sij, int L);

/*
 * Calculate the L2 norm of a matrix and do not forget to take the root ;)
 */
arma::mat calc_sij(arma::cube wij, int L);



/*
 * This returns the old score, that is L2 norm of wij matrix (minus APC if specified so)
 */
arma::mat calcHeuristicAPC(int L, arma::cube wij_3d, bool apc, int debug);

/*
 * Calculate the amino acid frequency of a given MSA.
 * Gaps are NOT ignored in order to normalize the counts.
 */
arma::vec compute_global_aa_freq(int N, int L, std::string msafilename, std::string seqweightsfilename, int verbose);


/*
 * Return all pairwise amino acid frequencies (counts) for one alignment
 * pseudo-counts and weighting are used
 * gaps are NOT ignored for normalization
 */
arma::cube compute_pairwise_freq(int N, int L, std::string msafilename, std::string seqweightsfilename, double pseudocount_admixture_regulator, bool return_aa_counts, int verbose);

/*
 * Compute amino acid frequencies for each column of an MSA
 *
 * Add pseudocounts from global amino acid frequencies in form of an admixture
 * that depends on the diversity (Neff) of the alignment:
 *
 * (1 - tau) * aa_freq_i + (tau * global_aa_freq)
 * with tau  =  pseudocount_admixture_regulator / (Neff + pseudocount_admixture_regulator)
 *
 * tau = 0 --> observed aa frequencies without pseudocounts
 *
 */
arma::mat compute_aa_freq_per_column(int N, int L, std::string msafilename, std::string seqweightsfilename, double pseudocount_admixture_regulator, bool return_aa_counts, int verbose);


/*
 * Compute the (weighted) number of non_gapped sequences per position
 */
arma::vec compute_Ni(int N, int L, std::string msafilename, std::string seqweightsfilename, int verbose);

/********************************************************************************************************************************
 *
 * Python-Boost Wrapper Functions
 *
 *********************************************************************************************************************************/
/*
 * Wrapper for python
 */
boost::python::list calcAPC_py(int L, boost::python::list matrix, int debug);

/*
 * Wrapper for python
 */
boost::python::list calcHeuristicAPC_py(int L, std::string brawfilename, bool apc, int debug);

/*
 * Return a list of all sequence weights (len N)
 */
boost::python::list calculate_seq_weights_py(int N, int L, std::string msafilename, double reweighting_threshold, int debug);


/*
 * Python Wrapper for read_msa
 */
boost::python::list read_msa_py(std::string msafilename, int N, int L , int debug);


/*
 * Return a list of sequence weights (len N)
 */
boost::python::list read_sequence_weights_py(std::string seqweightsfilename, int N, int debug);


/*
 * Python Wrapper for read_qij_msgpack_py
 */
boost::python::list read_qij_msgpack_py(std::string qijabfilename, int i, int j, int debug);

/*
 * Python Wrapper for compute_Ni
 */
boost::python::list compute_Ni_py(int N, int L, std::string msafilename, std::string seqweightsfilename, int verbose);


/*
 * Python Wrapper for compute_global_aa_freq
 */
boost::python::list compute_global_aa_freq_py(int N, int L, std::string msafilename, std::string seqweightsfilename, int verbose);


/*
 * Python Wrapper for compute_pairwise_freq
 */
boost::python::list compute_pairwise_freq_py(int N, int L, std::string msafilename, std::string seqweightsfilename, double pseudocount_admixture_regulator, bool return_aa_counts, int verbose);


/*
 * Python Wrapper for compute_aa_freq_per_column
 */
boost::python::list compute_aa_freq_per_column_py(int N, int L, std::string msafilename, std::string seqweightsfilename, double pseudocount_admixture_regulator, bool return_aa_counts, int verbose);

/*
 * Python Wrapper for read_Nij_msgpack_py
 */
boost::python::list read_Nij_msgpack_py(std::string qijabfilename, int L, int debug);


/*
 * Python Wrapper for calculating qijab
 */
void write_qijab_msgpack(	int N,
							int L,
							double lfactor,
							double lambda_w_fix,
							double pseudocount_admixture_regulator,
							std::string msafilename,
							std::string seqweightsfilename,
							std::string brawfilename,
							std::string qijabfilename,
							int verbose);



#endif /* CONTACTUTILS */

