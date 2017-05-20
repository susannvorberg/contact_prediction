//============================================================================
// Name        : boost_converters.cpp
// Author      : susi
// Version     :
// Copyright   : Your copyright notice
// Description : Functions converting std::vector constructs to
//				 boost::python::list equivalents that can be
//				 returned to python
//============================================================================


#include "boost_converters.hpp"
#include <boost/python.hpp>
#include <map>

/*
 * ------------------------------------------------------------------  From Boost to STD
 *
 */






/*
 * convert boost::python::dict to std::map<>
 */
std::map<std::string, std::vector<double> > fromPythonDict(boost::python::dict& py_dict){
	std::map<std::string, std::vector<double> > extracted_map;

	boost::python::list keys = py_dict.keys();
	for (int i = 0; i < len(keys); ++i) {
	   boost::python::extract<std::string> extracted_key(keys[i]);
	   if(!extracted_key.check()){
			std::cout<<"Key invalid, map might be incomplete"<<std::endl;
			continue;
	   }
	   std::string key = extracted_key;
	   boost::python::extract<boost::python::list > extracted_val(py_dict[key]);
	   if(!extracted_val.check()){
	   std::cout<<"Value invalid, map might be incomplete"<<std::endl;
			continue;
	   }
	   std::vector<double> vec = fromPythonList<double>(extracted_val);
	   extracted_map[key] = vec;
	}

	return extracted_map;
 }


/*
 * convert boost::python::dict to std::map<>
 */
std::map<std::string, std::string> fromPythonDictToString(boost::python::dict& py_dict){
	std::map<std::string, std::string > extracted_map;

	boost::python::list keys = py_dict.keys();
	for (int i = 0; i < len(keys); ++i) {
	   boost::python::extract<std::string> extracted_key(keys[i]);
	   if(!extracted_key.check()){
			std::cout<<"Key invalid, map might be incomplete"<<std::endl;
			continue;
	   }
	   std::string key = extracted_key;
	   boost::python::extract<std::string > extracted_val(py_dict[key]);
	   if(!extracted_val.check()){
	   std::cout<<"Value invalid, map might be incomplete"<<std::endl;
			continue;
	   }
	   std::string value = std::string(extracted_val);
	   extracted_map[key] = value;
	}

	return extracted_map;
 }

/*
 * convert boost::python::list to std::vector<>
 */
template <typename T>
std::vector<T> fromPythonList(boost::python::list l) {
	//get length of list
	boost::python::ssize_t len_list = boost::python::len(l);

	//initialize arma::vec
	std::vector<T> outvec(len_list);

	for(int i = 0; i < len_list; i++){
		outvec[i] = boost::python::extract<T>(l[i]);
	}

	return outvec;
}
template std::vector<double> fromPythonList<double>(boost::python::list l);
template std::vector<int> fromPythonList<int>(boost::python::list l);


/*
 * Convert a boost::python::list to arma::vec
 */
arma::vec  py_list_to_armavec(boost::python::list &l){

	//get length of list
	boost::python::ssize_t len_list = boost::python::len(l);

	//initialize arma::vec
	arma::vec vector(len_list, arma::fill::zeros);

	for(int i = 0; i < len_list; i++){
		vector(i) = boost::python::extract<double>(l[i]);
	}

	return vector;

}

/*
 * Convert a boost::python::list to std::vec of type std::string
 */
std::vector<std::string>  py_list_to_stdvecstring(boost::python::list &l){

	//get length of list
	boost::python::ssize_t len_list = boost::python::len(l);

	//initialize arma::vec
	std::vector<std::string> outvec(len_list);


	for(int i = 0; i < len_list; i++){
		outvec[i] = boost::python::extract<std::string>(l[i]);
	}

	return outvec;

}

std::vector<std::vector<double> > py_list_to_stdvecvec(boost::python::list &l){
	//get length of list 1
	boost::python::ssize_t len_list = boost::python::len(l);

	//get length of list 2
	boost::python::list list_element_zero = boost::python::extract<boost::python::list>(l[0]);
	boost::python::ssize_t len_list2 = boost::python::len(list_element_zero);

	//initialize arma::vec
	std::vector<std::vector<double> > std_vecvec(len_list, std::vector<double>(len_list2, 0));

	for(int i = 0; i < len_list; i++){
		boost::python::list list_element_i = boost::python::extract<boost::python::list>(l[i]);

		for(int j=0; j<len_list2; j++){
			std_vecvec[i][j] = boost::python::extract<double>(list_element_i[j]);
		}
	}

	return std_vecvec;

}

/*
 * Convert a boost::python::list to arma::mat
 */
arma::mat  py_list_to_armamat(boost::python::list &l){

	//get length of list 1
	boost::python::ssize_t len_list = boost::python::len(l);

	//get length of list 2
	boost::python::list list_element_zero = boost::python::extract<boost::python::list>(l[0]);
	boost::python::ssize_t len_list2 = boost::python::len(list_element_zero);


	//initialize arma::vec
	arma::mat matrix(len_list, len_list2, arma::fill::zeros);

	for(int i = 0; i < len_list; i++){
		boost::python::list list_element_i = boost::python::extract<boost::python::list>(l[i]);
		for(int j=0; j<len_list2; j++){
			matrix(i,j) = boost::python::extract<double>(list_element_i[j]);
		}
	}

	return matrix;

}


/*
 * ------------------------------------------------------------------  From  STD to Boost
 *
 */



/*
 * convert std::vector<> to boost::python::list
 */
template <typename T>
boost::python::list toPythonList(std::vector<T> vector) {

    typename std::vector<T>::iterator iter;
    boost::python::list list;
    for (iter = vector.begin(); iter != vector.end(); ++iter) {
        list.append(*iter);
    }
    return list;
}
template boost::python::list toPythonList(std::vector<double> vector);
template boost::python::list toPythonList(std::vector<int> vector);




/*
 * convert std::vec to boost::python::list
 * todo: use arma::vec iterator, then we do not need to convert to std
 */


boost::python::list std_string_vector_to_py_list(std::vector<std::string> &v){

	boost::python::list result;
	std::vector<std::string>::iterator it;
	for (it = v.begin(); it != v.end(); ++it){
	   result.append(*it);
	}
	return result;
}


/*
 * convert arma::mat to boost::python::list
 */
template <typename T>
boost::python::list std_vectorvector_to_py_list(std::vector<std::vector<T> > v){

	boost::python::list result;

    typename std::vector< std::vector<T> >::iterator it;
    typename std::vector<T>::iterator row_it;


	for (it = v.begin(); it != v.end(); ++it){

		boost::python::list row;
		for (row_it = it->begin(); row_it != it->end(); ++row_it){
			row.append(*row_it);
		}
		result.append(row);

	}
	return result;
}
template boost::python::list std_vectorvector_to_py_list(std::vector<std::vector<double> > vector);
template boost::python::list std_vectorvector_to_py_list(std::vector<std::vector<int> > vector);



/*
 * Convert arma::mat to boost::python::list
 * row-wise
 */
boost::python::list armamat_to_py_list(arma::mat &matrix){

	std::vector<std::vector<double> > matrix_py(matrix.n_rows);

	for (size_t i = 0; i < matrix.n_rows; ++i) {
		matrix_py[i] = arma::conv_to<std::vector<double> >::from(matrix.row(i));
	}
	return(std_vectorvector_to_py_list(matrix_py));

}

/*
 * Converts a C++ map to a python dict
 */
//template <class K, class V>
//boost::python::dict toPythonDict(std::map<K, V> map) {
//    typename std::map<K, V>::iterator iter;
//	boost::python::dict dictionary;
//	for (iter = map.begin(); iter != map.end(); ++iter) {
//		dictionary[iter->first] = iter->second;
//	}
//	return dictionary;
//}
boost::python::dict toPythonDict(std::map<std::string, std::vector<double> > map) {
	std::map<std::string, std::vector<double> >::iterator iter;
	boost::python::dict dictionary;
	for (iter = map.begin(); iter != map.end(); ++iter) {
		dictionary[iter->first] = iter->second;
	}
	return dictionary;
}



