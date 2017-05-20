//============================================================================
// Name        : contactutils.cpp
// Author      : susi
// Copyright   : Your copyright notice
// Description : Implementation file for
//				 Functions needed in python
//				 Made accessable through boost::python
//============================================================================


#include "contactutils.hpp"
#include "boost_converters.hpp"
#include "util_math.hpp"

#include <omp.h>

#include <ESBTL/default.h>
#include <ESBTL/selected_atom_iterator.h>
#include <ESBTL/atom_selectors.h>

#include <msgpack.hpp>
#include <armadillo>
#include <numeric>
#include <exception>
#include <iostream>
#include <sstream>
#include <fstream>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

/********************************************************************************************************************************
 *
 * General Fct
 *
 *********************************************************************************************************************************/



/*
 * Convert amino acid ASCII code to an index, with 0 = gap and 1 = A, 2 = C, etc.
 */
unsigned char aatoi(unsigned char aa) {
	if(!isalpha(aa)) { return 0; }

	aa = toupper(aa);
	if(aa < 65 || aa > 90) { return 0; }
	return AMINO_INDICES[aa - 65];
}

unsigned char  itoaa(unsigned char i) {
	return CHAR_INDICES[i];
}


/*
 * Given index i and index j of an LxL matrix
 * return the index ij of the vectorized upper triangular matrix
 * without diagonal !!!
 * row-wise !!!
 */
double get_linear_index(int i, int j, int n){

	return( (n*(n-1)/2) - (n-i)*((n-i)-1)/2 + j - i - 1);

}

/*
 * From a linear index, rowwise, no diagonal
 * reconstruct paired index
 */
arma::vec get_paired_index(int ij, int n){

	int i = n - 2 - floor(sqrt(-8*ij + 4*n*(n-1)-7)/2.0 - 0.5);
	int j = ij + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2;

	arma::vec indices(2);

	indices(0) = i;
	indices(1) = j;

	return indices;

}

/*
 * Get the size of a file
 */
long getFileSize(FILE *file)
{
	long lCurPos, lEndPos;
	lCurPos = ftell(file);
	fseek(file, 0, 2);
	lEndPos = ftell(file);
	fseek(file, lCurPos, 0);
	return lEndPos;
}


/*
 * Setup  a matrix of size n x n
 * and fill the upper triangle with values from input vector
 * values needs to be of size n*(n+1)/2
 */
arma::mat fill_matrix_from_triangular_indices(arma::vec values, int n, bool indices_upper){

	assert(values.n_elem == n*(n+1)/2);

	//in order to find non-zero indices 9either lower or upper)
	arma::mat m(n,n,arma::fill::ones);

	//find all upper(or lower) triangular elements that are non-zero
	arma::uvec ind(values.n_elem, arma::fill::zeros);
	if(indices_upper) 	ind = arma::find(arma::trimatu(m));
	else 				ind = arma::find(arma::trimatl(m));

	//set matrix elements back to zero
	m.fill(0);

	//insert values in upper(or lower) triangle
	m.elem(ind) = values;

	//symmatrize martrix
	m = arma::symmatu(m);

	return(m);
}



/********************************************************************************************************************************
 *
 * Gaussian Mixtures
 *
 *********************************************************************************************************************************/

/*
 * inefficient way of calculating gaussian density with covMatrix as input
*/
double mydmvnorm_cov_logspace(const arma::vec &x, const arma::vec &meanVector, const arma::mat &cov, bool logd = false){

	arma::mat precMatrix(400,400,arma::fill::zeros);

	//Precision Matrix by inverting CovMatrix
	try{
		precMatrix = arma::inv(cov);
	}catch (const std::exception& e){
		std::cout << "mydmvnorm_cov_logspace: Cov Matrix not invertible." << std::endl;
		return 0;
	}
	arma::vec diff = x-meanVector;

	std::cout << "precMatrix(0,0): " << precMatrix(0,0) << ", precMatrix(0,1): " << precMatrix(0,1) << ", precMatrix(0,2): " << precMatrix(0,2) << ", precMatrix(1,0): " << precMatrix(1,0) << ", precMatrix(1,1): " << precMatrix(1,1)<< std::endl;


	//log determinant
	arma::vec eigval;
	arma::mat eigvec;
	double log_det_precMatrix = 0;
	try{

		arma::eig_sym(eigval, eigvec, precMatrix);
		log_det_precMatrix = sum(log(eigval));

		if( ! std::isnormal(log_det_precMatrix)){
			std::cout << "computation of log determinant via sum(log(eigenvalues)) results in : " << log_det_precMatrix << std::endl;

			arma::mat L = arma::chol(precMatrix);
			log_det_precMatrix = 2 * arma::accu(arma::log(L.diag()));

			if( ! std::isnormal(log_det_precMatrix)){
				std::cout << "computation of log determinant via 2 * arma::accu(arma::log(L.diag())) results in" << log_det_precMatrix <<std::endl;
				return 0;
			}
		}
	}catch (const std::exception& e){
		std::cout << "mydmvnorm_cov_logspace: Not able to calculate log Determinant of Prec Matrix: eigenvalue decomposition failed." << std::endl;
		return 0;
	}


	std::cout << "log_det_precMatrix: " << log_det_precMatrix << std::endl;
	std::cout << "x.n_elem*log(2*PI): " << x.n_elem*log(2*PI) << std::endl;
	std::cout << "(diff.t() * arma::symmatu(precMatrix) * diff): " << (diff.t() * arma::symmatu(precMatrix) * diff) << std::endl;


	//0.5 * (log(|SIGMA^-1|) - d*log(2PI) - (x-mean) SIGMA^-1 (x-mean)
	double log_density = arma::as_scalar(0.5 * ( log_det_precMatrix - x.n_elem*log2pi - (diff.t() * arma::symmatu(precMatrix) * diff) ));

	std::cout << "log_density: " << log_density << std::endl;

	if (logd==false) {
		log_density=exp(log_density);
	}
	return(log_density);

}



/*
 * Fast Implementation of multivariate gaussian density after baysm R package
 * using Prec Matrix as input
 */

arma::vec mydmvnorm_prec_logspace_chol(const arma::mat &x,
											const arma::rowvec &mean,
											const arma::mat &prec,
											bool logd = false,
											int cores = 1) {
		//omp_set_num_threads(cores);
		int n = x.n_rows;
		int xdim = x.n_cols;
		arma::vec out(n);
		arma::mat rooti = arma::trimatu(arma::chol(prec)); //cholesky decomposition of covMatrix = L , inverse of L is inverse of cholesky decomposition of prec matrix
		double rootisum = arma::sum(log(rooti.diag())); // calculation of log determinant = 2 * sum(log(diag(chol(precmat))))
		double constants = -(xdim/2) * log2pi;
		#pragma omp parallel for schedule(static)
		for (int i=0; i < n; i++) {
			arma::vec z = rooti * arma::trans( x.row(i) - mean) ;
			out(i)      = constants - 0.5 * arma::sum(z%z) + rootisum;
		}

		if (logd==false) {
		out=exp(out);
		}
		return(out);
}

/*
 * Fast Implementation of multivariate gaussian density after baysm R package
 * using Cov Matrix as input
 */
arma::vec dmvnrm_arma_mc(const arma::mat &x,
													const arma::rowvec &mean,
													const arma::mat &sigma,
													bool logd = false,
													int cores = 1) {

    //omp_set_num_threads(cores);
    int n = x.n_rows;
    int xdim = x.n_cols;
    arma::vec out(n);
    arma::mat rooti(400,400,arma::fill::zeros);


    bool success = false;
    arma::mat R(400,400);
    success = arma::chol(R, arma::symmatu(sigma));

    /*
    double perturbance = 1.0;
    while(success == false){
    	perturbance += 1e-6;
    	std::cout << "dmvnrm_arma_mc: Cholesky Decomposition of CovMatrix failed. Perturb CovMatrix by " <<  perturbance << std::endl;
    	std::cout << "CovMatrix:  " <<  min(sigma.diag()) << " " << max(sigma.diag()) << std::endl;
    	success = arma::chol(R, (arma::symmatu(sigma) + arma::eye(sigma.n_rows,sigma.n_rows) * perturbance));
    }
	*/

    if(success){
    	rooti = arma::trans(arma::inv(arma::trimatu(R))); //cholesky decomposition of covMatrix = L , inverse of L is inverse of cholesky decomposition of prec matrix
		double rootisum = arma::sum(log(rooti.diag()));
		double constants = -(xdim/2) * log2pi;
		#pragma omp parallel for schedule(static)
		for (int i=0; i < n; i++) {
			arma::vec z = rooti * arma::trans( x.row(i) - mean) ;
			out(i)      = constants - 0.5 * arma::sum(z%z) + rootisum;
		}
    }else{
    	std::cout << "dmvnrm_arma_mc: Cholesky Decomposition of CovMatrix failed. Try armadillo inversion and eigenvalue decomposition..." << std::endl;
		#pragma omp parallel for
		for (int i=0; i < n; i++) {
			out(i)      = mydmvnorm_cov_logspace(x.row(i).t(), mean.t(), arma::symmatu(sigma));
		}
    }


    if (logd==false) {
        out=exp(out);
    }
    return(out);
}

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
arma::mat read_msa(std::string msafilename, int N, int L , int debug=0) {
	if(debug > 0 ) std::cout << "Read in Psicov Alignment... \n";

	//initialize variable
	arma::mat psicov(N,L);

	char buf[SEQ_BUFFER];

	FILE* msafile = fopen(msafilename.c_str(), "r");
	if (msafile == NULL) {
		throw std::runtime_error(msafilename.c_str());
	}
	else {
		for(int i=0; i < N; i++) {
			fgets(buf, SEQ_BUFFER, msafile);
			for(int j = 0; j < L; j++) {
				psicov(i,j) = aatoi( buf[j] );
			}
		}
		fclose(msafile);
	}



	return psicov;
}

/*
 * Read vector of sequence weights from file
 *
 * each element contains the weight of a sequence from the alignment
 */
arma::vec read_sequenceweights(std::string seqweightsfilename, int N, int debug=0){


	if(debug  > 0) std::cout << "Read in Seq Weights from " << seqweightsfilename << " ..." << std::endl;

	//initialize variable
	arma::vec seq_weights(N, arma::fill::zeros);

	char buf[SEQ_BUFFER];

	FILE* seqweightfile = fopen(seqweightsfilename.c_str(), "r");
	if (seqweightfile == NULL) {
		throw std::runtime_error(seqweightsfilename.c_str());
	}
	else {
		for(int i=0; i < N; i++) {
			fgets(buf, SEQ_BUFFER, seqweightfile);
			std::istringstream iss(buf);
			iss >> seq_weights(i);
		}
		fclose(seqweightfile);
	}

	if(debug  > 1) std::cout << "seq_weights(0): " << seq_weights(0) << ", seq_weights(1): " << seq_weights(1)  << ", seq_weights(2): " << seq_weights(2) << std::endl;


	return seq_weights;

}

/*
 * Read in qij from qijab file for specified i and j
 *
 * qijab and Nij are necessary to compute Hessian
 */
std::vector<double> read_qij_msgpack(std::string qijabfilename, int i, int j, int debug){

	if(debug  > 0) std::cout << "Read msgpack qij file..." << std::endl;

	std::ifstream file;
	std::stringstream in;

	try {

		file.exceptions(std::ios::failbit | std::ios::badbit);
		file.open(qijabfilename.c_str(), std::ios_base::in | std::ios_base::binary);

		boost::iostreams::filtering_istreambuf decompressor;
		decompressor.push(boost::iostreams::gzip_decompressor());
		decompressor.push(file);
		boost::iostreams::copy(decompressor, in);
	}
    catch(const boost::iostreams::gzip_error& e) {
    	std::cout << "Cannot open gzipped qijab file " + qijabfilename << ": " << e.what() << std::endl;
    	exit (EXIT_FAILURE);
    }

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

    //get qij for specified i and j
    std::stringstream ss;
    ss << "i:" << i << ",j:" << j;
    std::string id = ss.str();
    std::map<std::string, msgpack::object> ijpair_map;
    msgpack_map.at(id).convert(&ijpair_map);

    std::vector<double> qij(400);
    ijpair_map.at("qij").convert(&qij);

	return qij;
}

/*
 * Read in Nij from qijab file
 *
 * qijab and Nij are necessary to compute Hessian
 */
arma::mat read_Nij_msgpack(std::string qijabfilename, int L, int debug){


	if(debug  > 0) std::cout << "Read msgpack qij file..." << std::endl;

	//initialise output vectors, are not returned in this method, use pointers
	arma::cube qijab_3d(L,L,400, arma::fill::zeros);							//without diagonal as we do not need pairs (i,i)
	//arma::vec N_ij(L*(L-1)/2, arma::fill::zeros);			//without diagonal as we do not need pairs (i,i)
	arma::mat N_ij(L,L,arma::fill::zeros);

	std::ifstream file;
	std::stringstream in;

	try {

		file.exceptions(std::ios::failbit | std::ios::badbit);
		file.open(qijabfilename.c_str(), std::ios_base::in | std::ios_base::binary);

		boost::iostreams::filtering_istreambuf decompressor;
		decompressor.push(boost::iostreams::gzip_decompressor());
		decompressor.push(file);
		boost::iostreams::copy(decompressor, in);
	}
    catch(const boost::iostreams::gzip_error& e) {
    	std::cout << "Cannot open gzipped qijab file " + qijabfilename << ": " << e.what() << std::endl;
    	exit (EXIT_FAILURE);
    }

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

	for(int i=0; i<(L-1); i++){
		for(int j=i+1; j<L; j++){

			std::stringstream ss;
			ss << "i:" << i << ",j:" << j;
			std::string id = ss.str();
			std::map<std::string, msgpack::object> ijpair_map;
			msgpack_map.at(id).convert(&ijpair_map);
			if(debug  > 0) std::cout << "id: " << id <<  std::endl;

			double Nij;
			ijpair_map.at("N_ij").convert(&Nij);
			//N_ij(i) = Nij;
			N_ij(i,j) = Nij;

			//std::vector<double> qij(400);
			//ijpair_map.at("qij").convert(&qij);
			//qijab_3d.tube(i,j)=arma::vec(qij);
			//double sum_qijab_3d = arma::sum(arma::vec(qij));
			//double mean_qijab_3d = arma::mean(arma::vec(qij));


			if(debug  > 0) std::cout << "N_ij(i,j) " << N_ij(i,j) <<  std::endl;
			//if(debug  > 0) std::cout << "qijab_3d(i,j,0) " << qijab_3d(i,j,0)<<  ", sum: " << sum_qijab_3d << ", mean: " << mean_qijab_3d << std::endl;

		}
	}

	return N_ij;

}


/*
 * Read in braw file and extract all w_ij[1:20, 1:20]
 */
arma::cube  read_braw(std::string brawfilename, int L, int debug=0) {
	if(debug > 0 ) std::cout << "Read in braw File... \n";

	//initialize variables
	arma::cube w_ij3d(L, L, 400, arma::fill::zeros);


	//deserialize
	msgpack::unpacked msg;


	//find out whether braw file is zipped or not
	if (brawfilename.find("gz") != std::string::npos) {

		if(debug > 0 ) std::cout << "braw  compressed!" << '\n';

		std::ifstream file;
		std::stringstream in;
		try {

			file.exceptions(std::ios::failbit | std::ios::badbit);
			file.open(brawfilename.c_str(), std::ios_base::in | std::ios_base::binary);

			boost::iostreams::filtering_istreambuf decompressor;
			decompressor.push(boost::iostreams::gzip_decompressor());
			decompressor.push(file);
			boost::iostreams::copy(decompressor, in);
		}
		catch(const boost::iostreams::gzip_error& e) {
	    	std::cout << "Cannot unzip braw file " + brawfilename << ": " << e.what() << std::endl;
	    	exit (EXIT_FAILURE);
		}

		file.close();

		// Get the size of the file in bytes
		long fileSize = in.str().size();

		// Convert strstream to const char*
		const std::string& tmp = in.str();
		const char *buf = tmp.c_str();


		msgpack::unpack(&msg, buf, fileSize, 0);

	}else{

		if(debug > 0 ) std::cout << "braw is uncompressed!" << '\n';

		FILE *file = NULL;      // File pointer

		// Open the file in binary mode using the "rb" format string
		// This also checks if the file exists and/or can be opened for reading correctly
		if ((file = fopen(brawfilename.c_str(), "rb")) == NULL)  {
			std::cout << "Cannot open " + brawfilename << "\n";
			exit (EXIT_FAILURE);
		}


		// Get the size of the file in bytes
		long fileSize = getFileSize(file);

		// Allocate space in the buffer for the whole file
		char *buf = new char[fileSize];

		// Read the file in to the buffer
		fread(buf, fileSize, 1, file);
		fclose(file);

		msgpack::unpack(&msg, buf, fileSize, 0);

		delete [] buf;

	}



		//msgpack::object obj = msg.get();
		//std::cout << "obj type: " << obj.type << std::endl;
		//std::cout << "map type: " << msgpack::type::MAP << std::endl;
		//std::cout << "obj size: " << obj.via.map.size << std::endl;


		//convert as a map of msgpack objects
		std::map<std::string, msgpack::object> msgpack_map;
		msg.get().convert(&msgpack_map);

		//get msgpack object for "x_pair"
		// the x_pair msgpack object is of type MAP and has keys "i/j" with j>i
		// this amounts to ncol*(ncol-1)/2 values in this map
		msgpack::object xpair_obj = msgpack_map.at("x_pair");

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
			xpair_ij_map.at("x").convert(&x_vec);

			//set the wijab object
			arma::mat x_mat(arma::vec(x_vec).begin(),21,21);    //column-wise to matrix
			arma::mat x_mat_sub = x_mat.submat(0, 0, 19, 19);

			//shoudl yield ab in row-wise order (i = row; j = column)
			w_ij3d.tube(i,j) = arma::vectorise(x_mat_sub);    //column wise back to vector


		}

		return w_ij3d;
}



/*
 * Calculate the weights for each sequence based on pairwise amino acid comparisons
 */
arma::vec calculate_seq_weights(int N,int L, arma::mat psicov, double reweighting_threshold, int debug=0){

	if(debug > 0 ) std::cout << "Calculate sequence weights... " << std::endl;

	arma::vec seq_weights(N, arma::fill::zeros);

	uint64_t nij = N * (N + 1) / 2;

	//pairwise sequence comparisons
	#pragma omp parallel for
	for(uint64_t ij = 0; ij < nij; ij++) {

		// compute i and j from ij
		// http://stackoverflow.com/a/244550/1181102
		uint64_t i, j;									// we need uint64_t because of precision
		uint64_t ii = N * (N + 1) / 2 - 1 - ij;
		uint64_t K = floor((sqrt(8 * ii + 1) - 1) / 2);
		i = N - 1 - K;
		j = ij - N * i + i * (i + 1) / 2;

		//iterate over all columns
		int identical_positions = 0;
		int relative_L = 0;

		for(int k=0; k<L; k++){
			if(psicov(i,k) != 0 && psicov(j,k) != 0){
				identical_positions += (psicov(i,k) == psicov(j,k)) ? 1 : 0;
				relative_L++;
			}
		}

		//if sequences i and j have more than reweighting_threshold identical (non-gapped) columns: increase weight
		if((double)identical_positions/relative_L > reweighting_threshold){
			#pragma omp atomic
			seq_weights(i) += 1;
			#pragma omp atomic
			seq_weights(j) += 1;
		}

	}

	for(int i=0 ; i < N; i++){
		seq_weights(i) = 1.0 / (seq_weights(i) -1);
	}
	return seq_weights;
}


/*
 * Calculate the APC correction for a given matrix
 */
arma::mat calc_apc(arma::mat Sij, int L) {

	arma::mat apc(L, L, arma::fill::zeros);

	arma::rowvec row_mean 	= arma::mean(arma::symmatu(Sij));
	double meanSij 			= arma::mean(row_mean);

	#pragma omp parallel for
	for(int i=0; i < (L-1); i++){
		for(int j=(i+1); j < L; j++){
			apc(i,j) = row_mean(i) * row_mean(j) / meanSij;
		}
	}

	return(apc);
}

/*
 * Calculate the L2 norm of a matrix and do not forget to take the root ;)
 */
arma::mat calc_sij(arma::cube wij, int L){

	arma::mat S_ij(L,L,arma::fill::zeros);

	#pragma omp parallel for
	for(int i=0; i < (L-1); i++){
		for(int j=(i+1); j < L; j++){
			arma::vec row = wij.tube(i,j);
			S_ij(i,j) = arma::norm(row, 2);
		}
	}
    return(S_ij);
}



/*
 * This returns the old score, that is L2 norm of wij matrix (minus APC if specified so)
 */
arma::mat calcHeuristicAPC(int L, arma::cube wij_3d, bool apc, int debug=0){

	//calculate l2 norm
	arma::mat S_ij = calc_sij(wij_3d, L);

	//apply apc corection
	if(apc){
		if(debug > 0 ) std::cout << "Applying APC Correction" << std::endl;
		arma::mat apc_m = calc_apc(S_ij, L);
		S_ij -= apc_m;
	}

	return(S_ij);
}

/*
 * Calculate the amino acid frequency of a given MSA.
 * Gaps are NOT ignored in order to normalize the counts.
 */
arma::vec compute_global_aa_freq(int N, int L, std::string msafilename, std::string seqweightsfilename, int verbose=0){

	//read alignment
	arma::mat psicov = read_msa(msafilename,  N,  L , verbose);

	//read seqweight file
	arma::vec seq_weights = read_sequenceweights(seqweightsfilename,  N, verbose);

	// Note for all computations:
	// psicov elements range from 0(=gap),1,2...20(=aa)

	arma::vec freq(21, arma::fill::zeros);
	for(int l=0; l < L; l++){
		for(int n=0; n < N; n++){
			freq(psicov(n,l)) += seq_weights(n);	//psicov elements range from 0(=gap),1,2...20
		}
	}

	//remove gap counts
	if (verbose > 0) std::cout << "Percentage of gaps (ignored): " <<  freq(0) / arma::sum(freq)  << std::endl;
	//freq(0) = 0;

	//this gives us the weighted aa frequencies of the alignment
	freq /= arma::sum(freq);

	if (verbose > 0) std::cout << "Sum of global aa freq is: " <<  arma::sum(freq) << std::endl;

	return(freq);
}



/*
 * Return all pairwise amino acid frequencies (counts) for one alignment
 * pseudo-counts and weighting are used
 * gaps are NOT ignored for normalization
 */
arma::cube compute_pairwise_freq(int N, int L, std::string msafilename, std::string seqweightsfilename, double pseudocount_admixture_regulator, bool return_aa_counts, int verbose=0){
	//read alignment
	arma::mat psicov = read_msa(msafilename,  N,  L , verbose);

	//read seqweight file
	arma::vec seq_weights = read_sequenceweights(seqweightsfilename,  N, verbose);
	double Neff = arma::sum(seq_weights);

	//initialize pair_frequency_cube
	arma::cube pair_freq(L,L,400, arma::fill::zeros);
	arma::cube pair_counts(L,L,400, arma::fill::zeros);

	//admixture coefficient tau
	double tau = pseudocount_admixture_regulator / (Neff + pseudocount_admixture_regulator);
	double one_minus_tau_sq = (1-tau)*(1-tau);
	if (verbose > 0) std::cout << "pseudocount admixture coefficient tau: " << tau << "and (1-tau)^2: " << one_minus_tau_sq << std::endl;

	//-------------------------------------------------------------------------------------------
	//// Determine single amino acid frequencies with pseudo counts
	//// Normalization WITH gaps
	//-------------------------------------------------------------------------------------------

	arma::mat aa_freq_i(21, L, arma::fill::zeros);
	for(int n = 0; n < N; n++) {
		for(int i = 0; i < L; i++) {
			aa_freq_i(psicov(n,i), i) += seq_weights(n);
		}
	}

	//set counts for all gaps to 0
	//now we DO use gaps for normalisation
	//aa_freq_i.row(0).fill(0);

	//global aa frequencies,  normalized !!with!! gaps
	arma::vec global_aa_freq(21, arma::fill::zeros);
	global_aa_freq = arma::sum(aa_freq_i, 1);//sum all rows --> col vec
	global_aa_freq /= arma::accu(global_aa_freq);
	if (verbose > 0) std::cout << "Sum of global aa freq is: " <<  arma::sum(global_aa_freq) << std::endl;

	//column-wise aa frequencies, normalized !!with!! gaps
	arma::rowvec colsum = arma::sum(aa_freq_i, 0);//sum all cols --> row vec
	aa_freq_i.each_row() /= colsum;
	if (verbose > 0) std::cout << "Sum of aa freq for all columns is: " << arma::accu(aa_freq_i) << " == " << L << "(L)" << std::endl;



	//add pseudocounts
	arma::mat aa_freq_i_pseudocounts(21,L, arma::fill::zeros);
	for(int i=0; i<L; i++){
		//compute empirical frequencies with admixture of pseudocounts
		aa_freq_i_pseudocounts.col(i) = ((1 - tau) * aa_freq_i.col(i)) + (tau * global_aa_freq);
	}
	if (verbose > 0) std::cout << "Sum of aa freq with pseudocounts for all columns is: " << arma::accu(aa_freq_i_pseudocounts) << " == " << L << "(L)" << std::endl;



	//-------------------------------------------------------------------------------------------
	//// Determine pairwise amino acid frequencies with pseudo counts
	//-------------------------------------------------------------------------------------------
	if (verbose > 0) std::cout << "Start iterating over all pairs i<j" << std::endl;
	for(int i=0; i<(L-1); i++){
		for(int j=i+1; j<L; j++){

			arma::mat aa_freq_ij(21, 21, arma::fill::zeros);
			for(int n=0; n < N; n++){
				aa_freq_ij(psicov(n,i), psicov(n,j)) += seq_weights(n);
			}

			//set counts for all gaps to 0
			//now we DO use gaps for normalisation
			//aa_freq_ij.row(0).fill(0);
			//aa_freq_ij.col(0).fill(0);

			//vectorise column-wise - without pseudo counts
			arma::mat aa_freq_ij_sub = aa_freq_ij.submat( 1, 1, 20, 20); //exclude row and column #0: gaps
			pair_counts.tube(i,j) = arma::vectorise(aa_freq_ij_sub.t()); //row-wise  == col-wise of transposed matrix

			//normalize !!with!! gaps
			double N_ij = arma::accu(aa_freq_ij);
			aa_freq_ij /= N_ij;

			//pairwise frequencies with pseudo counts
			arma::mat aa_freq_ij_pseudocounts(21,21, arma::fill::zeros);
			for(int a=1; a < 21; a++){
				for(int b=1; b < 21; b++){
					aa_freq_ij_pseudocounts(a,b) = one_minus_tau_sq * (aa_freq_ij(a,b) - (aa_freq_i(a,i) * aa_freq_i(b,j)) ) + (aa_freq_i_pseudocounts(a,i) * aa_freq_i_pseudocounts(b,j));
				}
			}

			//vectorise column-wise - with pseudo counts
			aa_freq_ij_pseudocounts = aa_freq_ij_pseudocounts.submat( 1, 1, 20, 20); //exclude row and column #0: gaps
			pair_freq.tube(i,j) = arma::vectorise(aa_freq_ij_pseudocounts, 1); //row-wise

		}//end j
	}//end i

	if (return_aa_counts) {
		return pair_counts;
	}else{
		return pair_freq;
	}


}


/*
 * Compute amino acid frequencies for each column of an MSA (21 states including gap == 0! )
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
arma::mat compute_aa_freq_per_column(int N, int L, std::string msafilename, std::string seqweightsfilename, double pseudocount_admixture_regulator, bool return_aa_counts, int verbose=0){

	//read alignment
	arma::mat psicov = read_msa(msafilename,  N,  L , verbose);

	//read seqweight file
	arma::vec seq_weights = read_sequenceweights(seqweightsfilename,  N, verbose);
	double Neff = arma::sum(seq_weights);

	//admixture coefficient tau
	double tau = pseudocount_admixture_regulator / (Neff + pseudocount_admixture_regulator);

	// Note for all computations:
	// psicov elements range from 0(=gap),1,2...20(=aa)

	//-------------------------------------------------------------------------------------------------
	// determine global (weighted) aa frequencies, normalized !!!with!!! gaps
	// determine 	q0(x_i=a) --> empirical frequencies, normalized !!!with!!! gaps
	//-------------------------------------------------------------------------------------------------

	arma::mat aa_freq_i(21, L, arma::fill::zeros);
	for(int n = 0; n < N; n++) {
		for(int i = 0; i < L; i++) {
			aa_freq_i(psicov(n,i), i) += seq_weights(n);
		}
	}

	//set counts for all gaps to 0 --> all further normalization ignores gaps
	//now we DO use gaps for normalisation
	//aa_freq_i.row(0).fill(0);

	//global amino acid counts
	arma::vec global_aa_counts = arma::sum(aa_freq_i, 1);//sum all rows --> col vec

	//global aa frequencies,  normalized !!with!! gaps
	arma::vec global_aa_freq = global_aa_counts / arma::accu(global_aa_counts);

	if (verbose > 0) std::cout << "Sum of global aa freq is: " <<  arma::sum(global_aa_freq) << std::endl;

	//aa frequencies at each position, normalized !!with!! gaps
	aa_freq_i /= Neff;

	if (verbose > 0) std::cout << "Sum of aa freq for all columns is: " << arma::accu(aa_freq_i) << " == " << L << "(L)" << std::endl;


	//-------------------------------------------------------------------------------------------------
	//	determine q(x_i=a)  --> empirical frequencies, with admixture of pseudocounts
	//-------------------------------------------------------------------------------------------------

	arma::mat aa_freq_i_pseudocounts(21,L, arma::fill::zeros);
	for(int i=0; i<L; i++){
		//compute empirical frequencies with admixture of pseudocounts
		aa_freq_i_pseudocounts.col(i) = ((1 - tau) * aa_freq_i.col(i)) + (tau * global_aa_freq);
	}
	if (verbose > 0) std::cout << "Sum of aa freq with pseudocounts for all columns is: " << arma::accu(aa_freq_i_pseudocounts) << " == " << L << "(L)" << std::endl;

	//instead of the frequencies, return the amino acid counts (gaps ARE counted for freq calculation;)
	if(return_aa_counts){
		aa_freq_i_pseudocounts *= Neff;
	}

	return(aa_freq_i_pseudocounts);

}


/*
 * Compute the (weighted) number of non_gapped sequences per position
 */
arma::vec compute_Ni(int N, int L, std::string msafilename, std::string seqweightsfilename, int verbose=0){

	//read alignment
	arma::mat psicov = read_msa(msafilename,  N,  L , verbose);

	//read seqweight file
	arma::vec seq_weights = read_sequenceweights(seqweightsfilename,  N, verbose);
	double Neff = arma::sum(seq_weights);

	arma::vec Ni(L, arma::fill::zeros);
	for(int i = 0; i < L; i++) {

		arma::vec aa_freq_i(21, arma::fill::zeros);
		for(int n = 0; n < N; n++) {
			aa_freq_i(psicov(n,i)) += seq_weights(n);
		}
		Ni(i) = Neff - aa_freq_i(0);
	}


	return(Ni);

}



/********************************************************************************************************************************
 *
 * Python-Boost Wrapper Functions
 *
 *********************************************************************************************************************************/

/*
 * Wrapper for python
 */
boost::python::list calcAPC_py(int L, boost::python::list matrix, int debug=0){

	arma::mat m = py_list_to_armamat(matrix);

	arma::mat apc = calc_apc(m, L);

	return armamat_to_py_list(apc);

}


/*
 * Wrapper for python
 */
boost::python::list calcHeuristicAPC_py(int L, std::string brawfilename, bool apc, int debug=0){

	arma::cube wij_3d = read_braw(brawfilename, L, debug);

	arma::mat apc_score_matrix = calcHeuristicAPC(L, wij_3d, apc, debug);

	std::vector<std::vector<double> > apc_score_matrix_std(apc_score_matrix.n_rows);

	for (size_t i = 0; i < apc_score_matrix.n_rows; ++i) {
		apc_score_matrix_std[i] = arma::conv_to<std::vector<double> >::from(apc_score_matrix.row(i));
	}

	return(std_vectorvector_to_py_list(apc_score_matrix_std));

}


/*
 * Return a list of all sequence weights (len N)
 */
boost::python::list calculate_seq_weights_py(int N, int L, std::string msafilename, double reweighting_threshold, int debug=0){

	arma::mat psicov = read_msa(msafilename,  N,  L , debug);

	arma::vec seq_weights = calculate_seq_weights(N, L, psicov, reweighting_threshold, debug);

	std::vector<double> std_seq_weights = arma::conv_to<std::vector<double> >::from(seq_weights);

	return(toPythonList(std_seq_weights));
}

/*
 * Python Wrapper for read_msa
 */
boost::python::list read_msa_py(std::string msafilename, int N, int L , int debug=0){

	arma::mat msa = read_msa(msafilename, N, L , debug);

	//convert arma::mat to boost::python::list
	std::vector<std::vector<double> >msa_std(msa.n_rows);
	for (size_t i = 0; i < msa.n_rows; ++i) {
		msa_std[i] = arma::conv_to<std::vector<double> >::from(msa.row(i));
	}

	return(std_vectorvector_to_py_list(msa_std));
}



/*
 * Return a list of sequence weights (len N)
 */
boost::python::list read_sequence_weights_py(std::string seqweightsfilename, int N, int debug=0){

	arma::vec seq_weights = read_sequenceweights(seqweightsfilename, N, debug);
	std::vector<double> std_seq_weights = arma::conv_to<std::vector<double> >::from(seq_weights);
	return(toPythonList(std_seq_weights));
}

/*
 * Python Wrapper for compute_Ni
 */
boost::python::list compute_Ni_py(int N, int L, std::string msafilename, std::string seqweightsfilename, int verbose=0){

	arma::vec Ni = compute_Ni(N, L, msafilename, seqweightsfilename,  verbose);
	std::vector<double> std_Ni = arma::conv_to<std::vector<double> >::from(Ni);
	return(toPythonList(std_Ni));
}



/*
 * Python Wrapper for compute_global_aa_freq_py
 */
boost::python::list compute_global_aa_freq_py(int N, int L, std::string msafilename, std::string seqweightsfilename, int verbose=0){

	arma::vec global_aa_freq = compute_global_aa_freq(N, L, msafilename, seqweightsfilename, verbose);
	std::vector<double> std_global_aa_freq = arma::conv_to<std::vector<double> >::from(global_aa_freq);
	return(toPythonList(std_global_aa_freq));

}




/*
 * Python Wrapper for compute_pairwise_freq
 */
boost::python::list compute_pairwise_freq_py(int N, int L, std::string msafilename, std::string seqweightsfilename,  double pseudocount_admixture_regulator, bool return_aa_counts, int verbose){

	//cube with pairwise frequencies for all pairs i,j with j>i
	arma::cube pair_aa_freq = compute_pairwise_freq(N, L, msafilename, seqweightsfilename, pseudocount_admixture_regulator, return_aa_counts, verbose);

	if (verbose > 0) std::cout << "start converting cube to boost::python::list" << std::endl;


	boost::python::list result;

	for (int i=0; i<(L-1); i++ ){
		boost::python::list row_i;

		for (int j=i+1; j<L; j++){
			arma::vec pair_i_j = pair_aa_freq.tube(i,j);
			boost::python::list row_j;

			//ab is row-wise
			for (int x=0; x<400; x++){
				row_j.append(pair_i_j(x)); //third level
			}

			row_i.append(row_j); //second level
		}
		result.append(row_i); //first level
	}

	return result;
}


/*
 * Python Wrapper for compute_aa_freq_per_column
 */
boost::python::list compute_aa_freq_per_column_py(int N, int L, std::string msafilename, std::string seqweightsfilename, double pseudocount_admixture_regulator, bool return_aa_counts, int verbose=0){

	arma::mat aa_freq_per_column = compute_aa_freq_per_column(N, L, msafilename, seqweightsfilename, pseudocount_admixture_regulator,  return_aa_counts, verbose);

	//convert arma::mat to boost::python::list
	std::vector<std::vector<double> >aa_freq_per_column_std(aa_freq_per_column.n_rows);
	for (size_t i = 0; i < aa_freq_per_column.n_rows; ++i) {
		aa_freq_per_column_std[i] = arma::conv_to<std::vector<double> >::from(aa_freq_per_column.row(i));
	}

	return(std_vectorvector_to_py_list(aa_freq_per_column_std));

}

/*
 * Python Wrapper for read_qij_msgpack_py
 */
boost::python::list read_qij_msgpack_py(std::string qijabfilename, int i, int j, int debug){

	std::vector<double> qij = read_qij_msgpack(qijabfilename, i, j, debug);

	return(toPythonList(qij));

}



/*
 * Python Wrapper for read_Nij_msgpack_py
 */
boost::python::list read_Nij_msgpack_py(std::string qijabfilename, int L, int debug){

	arma::mat N_ij = read_Nij_msgpack(qijabfilename, L, debug);

	//convert arma::mat to boost::python::list
	std::vector<std::vector<double> >N_ij_std(N_ij.n_rows);
	for (size_t i = 0; i < N_ij.n_rows; ++i) {
		N_ij_std[i] = arma::conv_to<std::vector<double> >::from(N_ij.row(i));
	}

	return(std_vectorvector_to_py_list(N_ij_std));

}

/*
 * Calculating q'_ij and N_ij
 * see theory eq. 36
 */

void write_qijab_msgpack(int N,
		int L,
		double lfactor,
		double lambda_w_fix,
		double pseudocount_admixture_regulator,
		std::string msafilename,
		std::string seqweightsfilename,
		std::string brawfilename,
		std::string qijabfilename,
		int verbose = 0){

	//determine regularization constant
	double lambda_w			= lfactor * (L-1) + lambda_w_fix;

	//read braw
	arma::cube w_ij3d = read_braw(brawfilename, L, verbose);

	//read alignment
	arma::mat psicov = read_msa(msafilename,  N,  L , verbose);

	//read seqweight file
	arma::vec seq_weights = read_sequenceweights(seqweightsfilename,  N, verbose);
	double Neff = arma::sum(seq_weights);

	//admixture coefficient tau
	double tau = pseudocount_admixture_regulator / (Neff + pseudocount_admixture_regulator);
	double one_minus_tau_sq = (1-tau)*(1-tau);
	if (verbose > 0) std::cout << "pseudocount admixture coefficient tau: " << tau << "and (1-tau)^2: " << one_minus_tau_sq << std::endl;


	// Note for all computations:
	// psicov elements range from 0(=gap),1,2...20(=aa)

	//-------------------------------------------------------------------------------------------------
	// determine global (weighted) aa frequencies, normalized without gaps
	// determine 	q0(x_i=a) --> empirical frequencies, normalized without gaps
	//-------------------------------------------------------------------------------------------------
	arma::mat aa_freq_i(21, L, arma::fill::zeros);
	for(int n = 0; n < N; n++) {
		for(int i = 0; i < L; i++) {
			aa_freq_i(psicov(n,i), i) += seq_weights(n);
		}
	}

	//set counts for all gaps to 0
	aa_freq_i.row(0).fill(0);

	//global aa frequencies,  normalized without gaps
	arma::vec global_aa_freq(21, arma::fill::zeros);
	global_aa_freq = arma::sum(aa_freq_i, 1);//sum all rows --> col vec
	global_aa_freq /= arma::accu(global_aa_freq);

	if (verbose > 0) std::cout << "Sum of global aa freq is: " <<  arma::sum(global_aa_freq) << std::endl;

	//column-wise aa frequencies, normalized without gaps
	arma::rowvec colsum = arma::sum(aa_freq_i, 0);//sum all cols --> row vec
	aa_freq_i.each_row() /= colsum;

	if (verbose > 0) std::cout << "Sum of aa freq for all columns is: " << arma::accu(aa_freq_i) << " == " << L << "(L)" << std::endl;


	//-------------------------------------------------------------------------------------------------
	//	determine q(x_i=a)  --> empirical frequencies, with admixture of pseudocounts
	//-------------------------------------------------------------------------------------------------

	arma::mat aa_freq_i_pseudocounts(21,L, arma::fill::zeros);
	for(int i=0; i<L; i++){
		//compute empirical frequencies with admixture of pseudocounts
		aa_freq_i_pseudocounts.col(i) = ((1 - tau) * aa_freq_i.col(i)) + (tau * global_aa_freq);
	}
	if (verbose > 0) std::cout << "Sum of aa freq with pseudocounts for all columns is: " << arma::accu(aa_freq_i_pseudocounts) << " == " << L << "(L)" << std::endl;




	//Create output file object and buffer
	FILE *out = fopen(qijabfilename.c_str(), "w");
	msgpack::sbuffer buffer;
	msgpack::packer<msgpack::sbuffer> pk(&buffer);

	//pack L*(L+1)/2 new map objects
	pk.pack_map(L*(L-1)/2);


	for(int i=0; i<(L-1); i++){
		for(int j=i+1; j<L; j++){


			//-------------------------------------------------------------------------------------------------
			// determine q0(x_i=a, x_j=b) --> empirical frequencies, normalized without gaps
			// 			 q(x_i=a, x_j=b)  --> empirical frequencies, with pseudocounts
			//-------------------------------------------------------------------------------------------------
			arma::mat aa_freq_ij(21, 21, arma::fill::zeros);
			for(int n=0; n < N; n++){
				aa_freq_ij(psicov(n,i), psicov(n,j)) += seq_weights(n);
			}

			//set counts for all gaps to 0
			aa_freq_ij.row(0).fill(0);
			aa_freq_ij.col(0).fill(0);

			//normalize without gaps, eq 13 in theory
			double N_ij = arma::accu(aa_freq_ij);
			aa_freq_ij /= N_ij;

			arma::mat aa_freq_ij_pseudocounts(21,21, arma::fill::zeros);
			for(int a=1; a < 21; a++){
				for(int b=1; b < 21; b++){
					aa_freq_ij_pseudocounts(a,b) = one_minus_tau_sq * (aa_freq_ij(a,b) - (aa_freq_i(a,i) * aa_freq_i(b,j)) ) + (aa_freq_i_pseudocounts(a,i) * aa_freq_i_pseudocounts(b,j));
				}
			}

			aa_freq_ij_pseudocounts = aa_freq_ij_pseudocounts.submat( 1, 1, 20, 20); //exclude row and column #0: gaps
			arma::rowvec q_xia_xjb = arma::vectorise(aa_freq_ij_pseudocounts, 1); //dim=1 is row-wise

			//for debugging: should be 1
			if (verbose > 0) std::cout << "Sum of pairwise aa freq for columns " << i << " and " << j << ": " << arma::accu(aa_freq_ij) << std::endl;

			//for debugging: should be 1
			if (verbose > 0) std::cout << "Sum of pairwise aa freq with pseudocounts for columns " << i << " and " << j << ": " << arma::sum(q_xia_xjb) << std::endl;


			double lambda_w_div_Nij = lambda_w / N_ij;
			std::vector<double> vec_out(400);
			//this is row-wise ab! ! ! !
			for(int ab=0; ab < 400; ab++){
				vec_out[ab] 			= q_xia_xjb(ab) - lambda_w_div_Nij * w_ij3d(i,j,ab);
				if (vec_out[ab] < 0){
			    std::cout << "Warning: q'_ijab  is negative!" << std::endl;
			    std::cout << "i: " << i << ", j: " << j << ", ab:" << ab << ", vec_out[ab] : " << vec_out[ab]  << ", lambda_w: " << lambda_w << ", N_ij: " << N_ij << ", w_ij3d(i,j,ab): " << w_ij3d(i,j,ab) << ", q_xia_xjb(ab): " << q_xia_xjb(ab) << std::endl;
			    }
			}

/*          double min_vecout  = *std::min_element(vec_out.begin(), vec_out.end());
			if (min_vecout < 0){
			    std::cout << "Warning: q'_ij contains negative values!" << std::endl;
			    int pos = std::find(vec_out.begin(), vec_out.end(), min_vecout) - vec_out.begin();
			    std::cout << "i: " << i << ", j: " << j << ", ab:" << pos << ", vec_out[ab] : " << vec_out[pos]  << ", lambda_w: " << lambda_w << ", N_ij: " << N_ij << ", w_ij3d(i,j,ab): " << w_ij3d(i,j,pos) << ", q_xia_xjb(ab): " << q_xia_xjb(pos) << std::endl;
			}*/

			//this should be the case becaue sum_(a,b) q(x_i=a,x_j=b) = 1 and sum_(a,b) w_ij(a,b) = 0
			//double sum_vec_out = std::accumulate(vec_out.begin(), vec_out.end(), 0);
			//assert(sum_vec_out == 1 && "q'(ij) != 1 !")

			//pack map object with Nij and qij
			std::stringstream ss;
			ss << "i:" << i << ",j:" << j;
			std::string id = ss.str();
			pk.pack(id);
			pk.pack_map(2);

			pk.pack(std::string("N_ij"));
			pk.pack(N_ij);
			pk.pack(std::string("qij"));
			pk.pack(vec_out);


		}
	}

	fwrite(buffer.data(), buffer.size(), 1, out);
	fclose(out);

	//N_ij are not influenced by pseudo-counts; N_ij simply represents the number of non-gapped sequence pairs
	//the raw pair frequencies q(x_i=a, x_j=b) though, they are influenced by pseudo counts, such that q(x_i=a, x_j=b) == 0 will get a non_zero value

}


/********************************************************************************************************************************
 *
 * Definition of Boost Module
 *
 *********************************************************************************************************************************/



/*
 * Define the BOOST Module
 */
BOOST_PYTHON_MODULE(libcontactutils)
{

	boost::python::def("calcAPC_py", calcAPC_py);
	boost::python::def("calcHeuristicAPC_py", calcHeuristicAPC_py);
	boost::python::def("read_msa_py", read_msa_py);
	boost::python::def("calculate_seq_weights_py", calculate_seq_weights_py);
	boost::python::def("read_sequence_weights_py", read_sequence_weights_py);
	boost::python::def("write_qijab_msgpack", write_qijab_msgpack);
	boost::python::def("read_Nij_msgpack_py", read_Nij_msgpack_py);
	boost::python::def("read_qij_msgpack_py", read_qij_msgpack_py);
	boost::python::def("compute_global_aa_freq_py", compute_global_aa_freq_py);
	boost::python::def("compute_pairwise_freq_py", compute_pairwise_freq_py);
	boost::python::def("compute_aa_freq_per_column_py", compute_aa_freq_per_column_py);
	boost::python::def("compute_Ni_py", compute_Ni_py);

}
