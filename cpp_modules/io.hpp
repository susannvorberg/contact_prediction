//============================================================================
// Name        : io.hpp
// Author      : Susann Vorberg
// Copyright   : Your copyright notice
// Description : Header file for
//				 functions for file input/output
//				 Made accessable through boost::python
//============================================================================

#ifndef BAYESIAN_IO
#define BAYESIAN_IO


#include <boost/python.hpp>
#include <armadillo>


#define SEQ_BUFFER 8192
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
 * Get the size of a file
 */
long getFileSize(FILE *file);


/** Read an MSA file into an index matrix.
 *
 * Matrix will be returned as a pointer to a one-dimensional int array
 * where the cell at row i and column j will have the index i + nrow*j
 */
std::vector< std::vector<int> > read_msa(std::string msafilename , int debug);


void read_braw_meta(std::string brawfilename, int debug);


/*
 * Read in braw file and extract all w_ij[1:20, 1:20]
 */
arma::cube read_braw(std::string brawfilename, int L, int debug);



/*
 * Read in qij from qijab file
 *
 * qijab and Nij are necessary to compute Hessian
 */
arma::mat read_q_msgpack(std::string qijabfilename, int L, int debug);



/*
 * Read in q_ijab and Nij from qijab file
 *
 * qijab and Nij are necessary to compute Hessian
 */
arma::mat read_Nij_msgpack(std::string qijabfilename, int L, int debug);

/********************************************************************************************************************************
 *
 * Python-Boost Wrapper Functions
 *
 *********************************************************************************************************************************/



/*
 * Python Wrapper for read_msa
 */
boost::python::list read_msa_py(std::string msafilename, int debug);



/*
 * Return a list of sequence weights (len N)
 */
boost::python::list read_sequence_weights_py(std::string seq_weights_filename, int debug);

/*
 * Python Wrapper for read_qij_msgpack_py
 */
boost::python::list read_q_msgpack_py(std::string qijabfilename, int L, int debug);


/*
 * Python Wrapper for read_Nij_msgpack_py
 */
boost::python::list read_Nij_msgpack_py(const std::string& qijabfilename, int L, int debug);




#endif /* BAYESIAN_IO */