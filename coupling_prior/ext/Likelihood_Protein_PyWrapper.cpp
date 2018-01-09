//============================================================================
// Name        : Likelihood_Protein_PyWrapper.cpp
// Author      : Susann Vorberg
// Version     :
// Copyright   : Your copyright notice
// Description : python wrapper functions for Protein LL
//============================================================================

#include <boost/python.hpp>
#include "Likelihood_Protein.hpp"
#include "Parameters.hpp"
#include "boost_converters.hpp"
#include <armadillo>
#include <map>

using namespace boost::python;


/*
 * Wrapper Constructor
 */
Likelihood_Protein* Likelihood_Protein_PyWrapper(
    std::string protein_id,
    std::string brawfilename,
    std::string qijabfilename,
    boost::python::dict parameters,
    bool L_dependent
    ) {

	//save parameters as a std::map
	std::map<std::string, std::vector<double> > parameterMap = fromPythonDict(parameters);


	return new Likelihood_Protein(protein_id, brawfilename, qijabfilename, parameterMap, L_dependent);

}

/*
 * Return the responsibilities p(k | ij) for all components k for a specified pair i,j
 */
boost::python::list get_neg_log_likelihood_pairwise_py(const Likelihood_Protein& protein){


	arma::vec negLL_vec = protein.get_neg_log_likelihood_pairwise();
	std::vector<double> std_negLL_vec = arma::conv_to<std::vector<double> >::from(negLL_vec);
	return(toPythonList(std_negLL_vec));
}


/*
 * PYTHON Wrapper for setting the batch residue pairs
 */
void set_pairs_py(Likelihood_Protein& protein, boost::python::list &i_indices_in, boost::python::list &j_indices_in, boost::python::list &contacts_in){

	arma::uvec i_indices_vec = arma::conv_to< arma::uvec >::from(py_list_to_armavec(i_indices_in));
	arma::uvec j_indices_vec = arma::conv_to< arma::uvec >::from(py_list_to_armavec(j_indices_in));
	arma::uvec contacts_vec  = arma::conv_to< arma::uvec >::from(py_list_to_armavec(contacts_in));

	protein.set_pairs(i_indices_vec, j_indices_vec, contacts_vec);
}


/*
 * Python wrapper for C++ Class using boost_python
 *
 * name of the .so file must match boost module name
*/
BOOST_PYTHON_MODULE(libproteinll){

	boost::python::class_<Likelihood_Protein>("Likelihood_Protein")

    //add the Wrapper Constructor
    //calling Protein_LaplacianApprox_PyWrapper constructor results in Protein_LaplacianApprox object
    //because Protein_LaplacianApprox_PyWrapper calls Protein_LaplacianApprox constructor
    .def( "__init__", boost::python::make_constructor( &Likelihood_Protein_PyWrapper) )


    .def("set_pairs_py", set_pairs_py)
    .def("get_neg_log_likelihood_pairwise_py", get_neg_log_likelihood_pairwise_py)
	.def("compute_negLL", &Likelihood_Protein::compute_negLL)
    ;
}




