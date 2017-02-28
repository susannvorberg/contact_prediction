//============================================================================
// Name        : LaplaceanApproxRegularizer_PyWrapper.hpp
// Author      : susi
// Version     :
// Copyright   : Your copyright notice
// Description : header file for python wrapper functions for regularizer
//============================================================================

#ifndef REG_PYWRAP
#define REG_PYWRAP

#include <boost/python.hpp>
#include "Parameters.hpp"
#include "Regularizer.hpp"

/*
 * Wrapper Constructor
 */
Regularizer* Regularizer_PyWrapper(	boost::python::dict parameters_,
                                    double regularization_parameter_mu_,
                                    double regularization_parameter_diagonal_PrecMat_,
                                    int debug_
                                    );


/********************************************************************************************************************************
 *
 * PYTHON WRAPPER FOR Gradients
 *
 *********************************************************************************************************************************/

/*
 * Python Wrapper for "gradient_mu_comp_reg"
 *
 * Calculate the gradient of the regularizer for mu for a component c
 */
boost::python::list gradient_mu_comp_reg_py(Regularizer& regularizer, int c);


/*
 * Python Wrapper for "gradient_choleskyL_comp"
 *
 * Computes the gradient wrt the cholesky decomposed precision Matrix L
 * from the gradient of the precision matrix of either:
 * - precision matrix Lambda
 * - regularized precision matrix Lambda
 */
boost::python::list gradient_diagprecMat_comp_reg_py(Regularizer& regularizer, int c);




#endif // REG_PYWRAP

