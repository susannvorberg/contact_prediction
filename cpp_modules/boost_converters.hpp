//============================================================================
// Name        : boost_converters.hpp
// Author      : susi
// Version     :
// Copyright   : Your copyright notice
// Description : Header file for
//				 Functions converting std::vector constructs to
//				 boost::python::list equivalents that can be
//				 returned to python
//============================================================================

#ifndef BOOST_CONVERTERS
#define BOOST_CONVERTERS

#include <boost/python.hpp>
#include <armadillo>
#include <map>

/*
 * convert boost::python::dict to std::map<>
 */
std::map<std::string, std::vector<double> > fromPythonDict(boost::python::dict& py_dict);

/*
 * convert boost::python::dict to std::map<>
 */
std::map<std::string, std::string> fromPythonDictToString(boost::python::dict& py_dict);


/*
 * convert boost::python::list to std::vector<>
 */
template <typename T>
std::vector<T> fromPythonList(boost::python::list l);


/*
 * convert std::vector<> to boost::python::list
 */
template <typename T>
boost::python::list toPythonList(std::vector<T> vector);




/*
 * convert std::vec to boost::python::list
 * todo: use arma::vec iterator, then we do not need to convert to std
 */
boost::python::list std_string_vector_to_py_list(std::vector<std::string> &v);

/*
 * convert arma::mat to boost::python::list
 */
boost::python::list std_vectorvector_to_py_list(std::vector<std::vector<double> > &v);


/*
 * Convert a boost::python::list to arma::vec
 */
arma::vec  py_list_to_armavec(boost::python::list &l);

/*
 * Convert a boost::python::list to std::vec of type std::string
 */
std::vector<std::string>  py_list_to_stdvecstring(boost::python::list &l);

/*
 * Convert a boost::python::list to arma::mat
 */
arma::mat  py_list_to_armamat(boost::python::list &l);

std::vector<std::vector<double> > py_list_to_stdvecvec(boost::python::list &l);

/*
 * Convert arma::mat to boost::python::list
 * row-wise
 */
boost::python::list armamat_to_py_list(arma::mat &matrix);


/*
 * Converts a C++ map to a python dict
 */
boost::python::dict toPythonDict(std::map<std::string, std::vector<double> > map);




#endif /* BOOST_CONVERTERS*/
