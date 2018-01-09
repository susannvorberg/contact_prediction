//============================================================================
// Name        : io.cpp
// Author      : Susann Vorberg
// Copyright   : Your copyright notice
// Description : implementation file for
//				 functions for file input/output
//				 Made accessable through boost::python
//============================================================================


#include "io.hpp"
#include "boost_converters.hpp"

//#include <pybind11/pybind11.h>
//#include <pybind11/stl.h>

#include <boost/python.hpp>
#include <armadillo>
#include <algorithm>
#include <iterator> /* istream_iterator */
#include <msgpack.hpp>

#include <stdio.h>
#include <cstdio> /* ftell fseek */
#include <ctype.h>  /* isalpha toupper */
#include <exception>
#include <sstream>   /* stringstream */
#include <fstream>
#include <iostream> /* cout */


#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

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




/** Read an MSA file into an index matrix.
 *
 * Matrix will be returned as a pointer to a one-dimensional int array
 * where the cell at row i and column j will have the index i + nrow*j
 */
std::vector< std::vector<int> > read_msa(std::string msafilename, int debug=0) {

	if(debug > 0 ) std::cout << "Read in Psicov Alignment... \n";


    std::vector<std::string> psicovLines;
    std::ifstream infile(msafilename);

    if(!infile){
        std::cout << "Error opening Alignment file " <<  msafilename << std::endl;
        return(std::vector< std::vector<int> >());
    }

    std::copy(std::istream_iterator<std::string>(infile),
                std::istream_iterator<std::string>(),
                back_inserter(psicovLines));

    int N = psicovLines.size();
    int L = psicovLines[0].size();


    if(debug > 0 ) std::cout << "Alignment has " << N << " sequences of length " << L << std::endl;

    //initialize psicov  vector of vectors
    std::vector< std::vector<int> > psicov(N, std::vector<int>(L));

    for(int n = 0; n < N; n++) {
        for(int l = 0; l < L; l++) {
            psicov[n][l] = aatoi( psicovLines[n][l] );
        }
    }

	return psicov;
}


//import bayesian_framework.cpp_functions.build.libio as io
//io.read_braw_meta('/home/vorberg/work/data//benchmarkset_cathV4/benchmarkset_cathV4_combs/ccmpred_dev_center_v/l_1772/braw/16vp_A_00.braw.gz', 5, 2)
void read_braw_meta(std::string brawfilename, int debug=0){
    if(debug > 0 ) std::cout << "Read in braw File... \n";


	//deserialize
	msgpack::unpacked msg;

	//find out whether braw file is zipped or not
	if (brawfilename.find("gz") != std::string::npos) {

		if(debug > 0 ) std::cout << "braw  compressed!" << '\n';

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

		if(debug > 0 ) std::cout << "braw is uncompressed!" << '\n';
		exit (EXIT_FAILURE);
    }


		//convert as a map of msgpack objects
		std::map<std::string, msgpack::object> msgpack_map;
		msg.get().convert(&msgpack_map);

    	//braw.meta['workflow'][0]['parameters']['regularization']['lambda_pair']
        std::map<std::string, msgpack::object> meta;
        msgpack_map.at("meta").convert(&meta);
        msgpack::object workflow = meta.at("workflow");
        std::vector<msgpack::object> workflow_list;
        workflow.convert(&workflow_list);
        std::map<std::string, msgpack::object> first_element;
        workflow_list[0].convert(&first_element);
        std::map<std::string, msgpack::object> parameters;
        first_element.at("parameters").convert(&parameters);
        std::map<std::string, msgpack::object> regularization;
        parameters.at("regularization").convert(&regularization);
        double lambda_pair;
        regularization.at("lambda_pair").convert(&lambda_pair);
        std::cout << "lambda_pair: " << lambda_pair << std::endl;

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

		if(debug > 0 ) std::cout << "braw is uncompressed!" << '\n';
		exit (EXIT_FAILURE);
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
			arma::mat x_mat(arma::vec(x_vec).begin(),21,21);    //length = 21, column-wise to matrix
			arma::mat x_mat_sub = x_mat.submat(0, 0, 19, 19);   //first column 0, last column 19

			//shoudl yield ab in row-wise order (i = row; j = column)
			w_ij3d.tube(i,j) = arma::vectorise(x_mat_sub);    //column wise back to vector


		}

		return w_ij3d;
}



/*
 * Read vector of sequence weights from file
 *
 * each element contains the weight of a sequence from the alignment
 */
//std::vector<double> read_sequence_weights(std::string seq_weights_filename, int debug=0)
//{
//
//
//	if(debug  > 0) std::cout << "Read in Seq Weights from " << seq_weights_filename << " ..." << std::endl;
//
//    std::vector<double> std_seq_weights;
//    std::ifstream ifs(seq_weights_filename);
//    if (ifs) {
//        std::copy(std::istream_iterator<double>(ifs),
//                    std::istream_iterator<double>(),
//                    std::back_inserter(std_seq_weights));
//    }
//    else {
//        std::cerr << "Couldn't open " << seq_weights_filename << " for reading\n";
//    }
//
//
//	if(debug  > 1) std::cout << "seq_weights(0): " << std_seq_weights[0] << ", seq_weights(1): " << std_seq_weights[1]  << ", seq_weights(2): " << std_seq_weights[2] << std::endl;
//
//	return std_seq_weights;
//}


/*
 * Read in qij from qijab file
 *
 * qijab and Nij are necessary to compute Hessian
 */
arma::mat read_q_msgpack(std::string qijabfilename, int L, int debug=0)
{
    if(debug  > 0) std::cout << "Read msgpack qij file..." << std::endl;

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
    	std::cout << "Cannot open gzipped qijab file " + qijabfilename << ": " << e.what() << std::endl;
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

	//read map object with keys N_ij and q_ij
	std::map<std::string, msgpack::object> msgpack_map;
	msg.get().convert(&msgpack_map);

    //read vector of size L*(L-1)/2
    msgpack::object q_ij_array = msgpack_map.at("q_ij");
    std::vector<double> qij_vec_std;
    q_ij_array.convert(&qij_vec_std);

    arma::mat qij_mat_Copy(&qij_vec_std.front(), 400, L*(L-1)/2);

    return qij_mat_Copy;


}


/*
 * Read in Nij from qijab file
 *
 * qijab and Nij are necessary to compute Hessian
 */
arma::mat read_Nij_msgpack(std::string qijabfilename, int L, int debug)
{


	if(debug  > 0) std::cout << "Read msgpack qij file..." << std::endl;

	//initialise output vectors, are not returned in this method, use pointers
	arma::cube qijab_3d(L,L,400, arma::fill::zeros);							//without diagonal as we do not need pairs (i,i)

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
    	std::cout << "Cannot open gzipped qijab file " + qijabfilename << ": " << e.what() << std::endl;
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

    //access of array at map key "N_ij"
	std::vector<double> Nij_std;
    msgpack_map.at("N_ij").convert(&Nij_std);
    arma::vec Nij_vec = arma::conv_to< arma::vec >::from(Nij_std);

    //fill upper triangle of matrix (without diagonal) with values
    arma::mat N_ij(L,L,arma::fill::ones);
    N_ij.diag().zeros();
    arma::uvec indices = arma::find(arma::trimatu(N_ij) > 0);
    N_ij.elem(indices) = Nij_vec;


	return N_ij;

}

/********************************************************************************************************************************
 *
 * Python-Boost Wrapper Functions
 *
 *********************************************************************************************************************************/





/*
 * Python Wrapper for read_msa
 */
boost::python::list read_msa_py(std::string msafilename , int debug=0)
{

	std::vector< std::vector<int> > msa = read_msa(msafilename , debug);

	return(std_vectorvector_to_py_list(msa));
}



/*
 * Python Wrapper for read_qij_msgpack_py
 */
boost::python::list read_q_msgpack_py(std::string qijabfilename, int L, int debug)
{

	arma::mat q = read_q_msgpack(qijabfilename, L, debug);

	return(armamat_to_py_list(q));

}


/*
 * Python Wrapper for read_Nij_msgpack_py
 */
boost::python::list read_Nij_msgpack_py(const std::string& qijabfilename, int L, int debug)
{

	arma::mat N_ij = read_Nij_msgpack(qijabfilename, L, debug);

    return(armamat_to_py_list(N_ij));

}


/********************************************************************************************************************************
 *
 * Definition of Boost Module
 *
 *********************************************************************************************************************************/


/*
 * Define the BOOST Module
 */
BOOST_PYTHON_MODULE(libio)
{
    boost::python::def("read_braw_meta", &read_braw_meta);
	boost::python::def("read_msa_py", &read_msa_py);
	boost::python::def("read_q_msgpack_py", &read_q_msgpack_py);
	boost::python::def("read_Nij_msgpack_py", &read_Nij_msgpack_py);

}

//namespace py = pybind11;
//
//PYBIND11_PLUGIN(libio) {
//    py::module m("libio", "library for reading and writing");
//
//    m.def("read_msa", &read_msa, "A function which reads a multiple sequence alignment",
//    py::arg("msafilename") = "", py::arg("debug") = 0);
//
//    m.def("read_Nij", &read_Nij_py, "A function which reads a Nij from a msgpacked qijab file");
//
//    return m.ptr();
//}
