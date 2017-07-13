//============================================================================
// Name        : Regularizer_PyWrapper.cpp
// Author      : susi
// Version     :
// Copyright   : Your copyright notice
// Description : python wrapper functions for regularizer of LL
//============================================================================

#include <boost/python.hpp>
#include "Regularizer.hpp"
#include "Regularizer_PyWrapper.hpp"
#include "boost_converters.hpp"
#include "Parameters.hpp"
#include <armadillo>
#include <map>

using namespace boost::python;


/*
 * Wrapper Constructor
 */
Regularizer* Regularizer_PyWrapper(
        boost::python::dict parameters_,
        double regularization_parameter_mu_,
        double regularization_parameter_diagonal_PrecMat_,
        bool L_dependent
    )
{

	//convert parameters from Py dict to std::map
	std::map<std::string, std::vector<double> > parameterMap = fromPythonDict(parameters_);

	return new Regularizer( parameterMap,
	                        regularization_parameter_mu_,
	                        regularization_parameter_diagonal_PrecMat_,
	                        L_dependent);

}



/*
 * Python wrapper for C++ Class using boost_python
 *
 * name of the .so file must match boost module name
*/
BOOST_PYTHON_MODULE(libreg){

	boost::python::class_<Regularizer>("Regularizer")

    //add the Wrapper Constructor
    .def( "__init__", boost::python::make_constructor( &Regularizer_PyWrapper) )

    .def("set_regularization_parameters", &Regularizer::set_regularization_parameters)
    .def("print_regularization_parameters", &Regularizer::print_regularization_parameters)

    //regularizers
    .def("regularizer_mu", &Regularizer::regularizer_mu)
    .def("regularizer_diagPrecMat", &Regularizer::regularizer_diagPrecMat)

    //gradients
    .def("gradient_diagprecMat_comp_reg", gradient_diagprecMat_comp_reg_py)
    .def("gradient_mu_comp_reg", gradient_mu_comp_reg_py)
    ;
}



/********************************************************************************************************************************
 *
 * PYTHON WRAPPER FOR GRADIENT
 *
 *********************************************************************************************************************************/


/*
 * Python Wrapper for "gradient_mu_comp_reg"
 *
 * Calculate the gradient of the regularizer for mu for a component c
 */
boost::python::list gradient_mu_comp_reg_py(Regularizer& regularizer, int c){

	std::vector<double> grad_mu_c = arma::conv_to<std::vector<double> >::from(regularizer.gradient_mu_comp_reg(c));
	return(toPythonList(grad_mu_c));
}



/*
 * Python Wrapper for "gradient_choleskyL_comp"
 *
 * Computes the gradient wrt the cholesky decomposed precision Matrix L
 * from the gradient of the precision matrix of either:
 * - precision matrix Lambda
 * - regularized precision matrix Lambda
 */
boost::python::list gradient_diagprecMat_comp_reg_py(Regularizer& regularizer, int c){

	arma::mat grad_diagprecMat_c_arma = regularizer.gradient_diagprecMat_comp_reg(c);
	return(armamat_to_py_list(grad_diagprecMat_c_arma));
}







