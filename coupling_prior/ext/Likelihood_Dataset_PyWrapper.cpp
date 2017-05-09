//============================================================================
// Name        : Likelihood_Dataset_PyWrapper.cpp
// Author      : Susann Vorberg
// Version     :
// Copyright   : Your copyright notice
// Description : python wrapper functions for LL
//============================================================================

#include <boost/python.hpp>
#include "Likelihood_Dataset.hpp"
#include "Likelihood_Protein.hpp"
#include "Parameters.hpp"
#include "boost_converters.hpp"
#include <armadillo>
#include <map>

using namespace boost::python;



/*
 * Python wrapper for C++ Class using boost_python
 *
 * name of the .so file must match boost module name
*/
BOOST_PYTHON_MODULE(libll){

	boost::python::class_<Likelihood_Dataset>("Likelihood_Dataset",
	                                            boost::python::init<boost::python::dict, boost::python::dict>())

    .def("get_f", &Likelihood_Dataset::get_f)
	.def("get_gradient_dict", &Likelihood_Dataset::get_gradient_dict)

	.def("set_debug_mode", &Likelihood_Dataset::set_debug_mode)
	.def("set_threads_per_protein", &Likelihood_Dataset::set_threads_per_protein)
	.def("set_parameters", &Likelihood_Dataset::set_parameters)


	.def("compute_f_df", &Likelihood_Dataset::compute_f_df)
	.def("compute_f", &Likelihood_Dataset::compute_f)
    ;
}







