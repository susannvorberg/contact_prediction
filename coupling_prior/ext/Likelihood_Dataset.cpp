//============================================================================
// Name        : Likelihood_Dataset.cpp
// Author      : susi
// Version     :
// Copyright   : Your copyright notice
// Description : return neg ll and gradients for dataset of proteins
//============================================================================


#include "Likelihood_Dataset.hpp"
#include "Likelihood_Protein.hpp"
#include "boost_converters.hpp"
#include <boost/python.hpp>
#include <omp.h>
#include <armadillo>
#include <map>


/*
 * Constructor
 */

Likelihood_Dataset::Likelihood_Dataset(
                    boost::python::dict pairs_per_protein_,
                    boost::python::dict parameters_
                    )
{

    //transform dict to map of protein structs
	protein_dict_to_protein_map(pairs_per_protein_);
	nr_proteins = this->dataset.size();

	debug_mode		= 0;
    nr_threads_prot = 1;

	//save parameters as a std::map
	parameterMap = fromPythonDict(parameters_);
	nr_components = parameterMap.size()/ 4 ;



	//initialize variables
	f = 0.0;
	grad_weight_bg.zeros(nr_components);
	grad_weight_contact.zeros(nr_components);
	grad_mu.zeros(400, nr_components);
    grad_precMat.zeros(400,nr_components);

    if (debug_mode > 0){
        std::cout << "Initialised LL Dataset class with following settings: " << std::endl;
        std::cout << "debug_mode: " << debug_mode << ", nr_components:" << nr_components << ", nr_threads:" << nr_threads_prot<< std::endl;
    }

}



/*
 * Transform the dictionary/list from python
 * into a map of std::vec of protein structs
 */
void Likelihood_Dataset::protein_dict_to_protein_map(boost::python::dict &pairs_per_protein)
{

	boost::python::list keys = pairs_per_protein.keys();
	dataset = std::vector<myProtein>(len(keys));
	for (int i = 0; i < len(keys); i++) {

		boost::python::dict protein_dict = boost::python::extract<boost::python::dict>(pairs_per_protein[keys[i]]);

		//protein stats
		std::string protein_name = boost::python::extract<std::string>(keys[i]);
		int N = boost::python::extract<int>(protein_dict["N"]);
		int L = boost::python::extract<int>(protein_dict["L"]);
		std::string braw_file_path  = boost::python::extract<std::string>(protein_dict["braw_file_path"]);
		std::string qijabfilename   = boost::python::extract<std::string>(protein_dict["qijabfilename"]);

		boost::python::list residue_i_list = boost::python::extract<boost::python::list>(protein_dict["residue_i"]);
		boost::python::list residue_j_list = boost::python::extract<boost::python::list>(protein_dict["residue_j"]);
		boost::python::list contact_list   = boost::python::extract<boost::python::list>(protein_dict["contact"]);

	    arma::uvec residue_i    = arma::conv_to< arma::uvec >::from(py_list_to_armavec(residue_i_list));
	    arma::uvec residue_j    = arma::conv_to< arma::uvec >::from(py_list_to_armavec(residue_j_list));
	    arma::uvec contact      = arma::conv_to< arma::uvec >::from(py_list_to_armavec(contact_list));

		myProtein protein {
		   protein_name,
		   N,
		   L,
		   braw_file_path,
		   qijabfilename,
		   residue_i,
		   residue_j,
		   contact
		};

		dataset[i] = protein;

	}

}


/*
* Set the level of printing debugging information
*/
void Likelihood_Dataset::set_debug_mode(int debug_mode){
    this->debug_mode = debug_mode;
}

/*
* Set the number of threads for parallelization
*/
void Likelihood_Dataset::set_threads_per_protein(int threads){
    this->nr_threads_prot = threads;
}


/*
* Set the parameters
*/
void Likelihood_Dataset::set_parameters(boost::python::dict parameters){
	this->parameterMap = fromPythonDict(parameters);
	this->nr_components = parameterMap.size()/ 4 ;
}


/*
 * return the neg LL of the protein dataset according to specified parameters
 */
double Likelihood_Dataset::get_f()
{
	return(this->f);
}


/*
 * return the gradients for all parameters in the likelihood function
 * according to order of parameter names list
 */
boost::python::dict Likelihood_Dataset::get_gradient_dict()
{

	boost::python::dict gradient_dict;

	typedef std::map<std::string, std::vector<double> >::iterator it_type;
	for(it_type iterator = parameterMap.begin(); iterator != parameterMap.end(); iterator++) {

	    std::string parameter_name = iterator->first;

	    //determine component number
		int position = parameter_name.rfind ('_');
		int component = std::stoi(parameter_name.substr(position+1));

	    if (parameter_name.find("weight_bg") != std::string::npos) {
			boost::python::list grad;
			grad.append(grad_weight_bg(component));

			gradient_dict[parameter_name] = grad;
		}
		else if (parameter_name.find("weight_contact") != std::string::npos) {
			boost::python::list grad;
			grad.append(grad_weight_contact(component));

			gradient_dict[parameter_name] = grad;
		}
		else if(parameter_name.find("mu") != std::string::npos){
			std::vector<double> std_grad_mu = arma::conv_to<std::vector<double> >::from(grad_mu.col(component));
			boost::python::list grad = toPythonList(std_grad_mu);

			gradient_dict[parameter_name] = grad;

		}
		else if(parameter_name.find("prec") != std::string::npos){

		    std::vector<double> std_grad_precMat = arma::conv_to<std::vector<double> >::from(grad_precMat.col(component));

			boost::python::list grad = toPythonList(std_grad_precMat);

			gradient_dict[parameter_name] = grad;
		}
	}

	return gradient_dict;
}

/*
 *  Actual computation of neg loglikelihood
 *  NO computation of gradients
 */
void Likelihood_Dataset::compute_f(){

	arma::vec f_vec(nr_proteins, arma::fill::zeros);


	if (this->debug_mode > 0){
	    std::cout << "Settings for OMP parallelization: " << std::endl;
	    std::cout << "#threads: " << nr_threads_prot  << " is_nested: " << omp_get_nested( ) << " is_dynamic: " << omp_get_dynamic( )<< " nr available cpus: " << omp_get_num_procs()<< std::endl;
	}


	//iterate over protein list
	#pragma omp parallel for schedule(dynamic) num_threads(nr_threads_prot)
	for(int p =0; p < nr_proteins; p++){

		myProtein protein_data = dataset[p];

		//create protein object
		if (this->debug_mode > 0){
    		#pragma omp critical
	        std::cout << "Compute protein " <<  protein_data.name << " (thread nr "<< omp_get_thread_num() << " on cpu " << sched_getcpu()  << ")" << std::endl;
        }

		Likelihood_Protein protein(
		        protein_data.name,
		        protein_data.N,
                protein_data.L,
                protein_data.brawfilename,
                protein_data.qijabfilename,
                this->parameterMap
                );

        //set pair information
		protein.set_pairs(protein_data.residue_i, protein_data.residue_j, protein_data.contact);


		//calc Hessian, mu_ij_k, lambda_ij_k, responsibilities
		//AND likelihood
		//with 1 thread
		protein.compute_negLL(1);

		//save protein likelihood value
		f_vec(p) = protein.get_neg_log_likelihood();

    }//end proteins

	//sum over all proteins (parallelize over proteins)
	f = arma::accu(f_vec);
}


/*
 * Actual computation of neg loglikelihood and gradients for whole dataset
 */
void Likelihood_Dataset::compute_f_df(int hessian_pseudocount)
{


    arma::vec f_vec(nr_proteins, arma::fill::zeros);

	arma::mat grad_weight_contact_mat(nr_proteins, nr_components, arma::fill::zeros);
	arma::mat grad_weight_bg_mat(nr_proteins, nr_components, arma::fill::zeros);
	arma::cube grad_mu_cube(nr_proteins, nr_components,400, arma::fill::zeros);
	arma::cube grad_precMat_cube(nr_proteins, nr_components,400, arma::fill::zeros);

	if (this->debug_mode > 0){
	    std::cout << "Settings for OMP parallelization: " << std::endl;
	    std::cout << "#threads: " << nr_threads_prot  << " is_nested: " << omp_get_nested( ) << " is_dynamic: " << omp_get_dynamic( )<< " nr available cpus: " << omp_get_num_procs() << "\n" << std::endl;
	}

    //iterate over protein list
	#pragma omp parallel for schedule(dynamic) num_threads(nr_threads_prot)
	for(int p =0; p < nr_proteins; p++){

		myProtein protein_data = dataset[p];

		//create protein object
		#pragma omp critical
		std::cout << "Compute protein " <<  protein_data.name << " (thread nr "<< omp_get_thread_num() << " on cpu " << sched_getcpu()  << ")" << std::endl;


		Likelihood_Protein protein(
		        protein_data.name,
		        protein_data.N,
                protein_data.L,
                protein_data.brawfilename,
                protein_data.qijabfilename,
                this->parameterMap
                );

        //set pair information
		protein.set_pairs(protein_data.residue_i, protein_data.residue_j, protein_data.contact);


		//calc Hessian, mu_ij_k, lambda_ij_k, responsibilities
		//likelihood and gradients for all pairs
		protein.compute_f_df(hessian_pseudocount);


		//save protein likelihood value
		f_vec(p) = protein.get_neg_log_likelihood();



        // get ALL gradients for ALL components
        for(int component = 0 ; component < this->nr_components; component++){

		    //weights
			grad_weight_contact_mat(p,component)    = protein.get_gradient_weight_comp(component, 1);
			grad_weight_bg_mat(p,component)         = protein.get_gradient_weight_comp(component, 0);
			//mu
			grad_mu_cube.tube(p,component)          = protein.get_gradient_mu_comp(component);
			//precMat - diagonal or isotrope
			grad_precMat_cube.tube(p,component)     = protein.get_gradient_precisionMatrix_comp(component);
            //precMat - isotrope (dependentL)
            //grad_precMat_cube.tube(p,component)     = protein.get_gradient_precisionMatrix_isotrop_Ldependent(component);
		}

	}//end proteins

	//sum over all proteins (parallelize over proteins)
	f = arma::accu(f_vec);

	arma::rowvec sum_grad_weight_contact    = arma::sum(grad_weight_contact_mat, 0); //sum of elements in each column (dim=0)
	arma::rowvec sum_grad_weight_bg         = arma::sum(grad_weight_bg_mat, 0);      //sum of elements in each column (dim=0)
	grad_weight_contact = arma::vec(sum_grad_weight_contact.t());
	grad_weight_bg		= arma::vec(sum_grad_weight_bg.t());

	arma::mat mat_grad_m 		= arma::sum(grad_mu_cube, 0);       //dim=0 indicates the sum of elements in each column within each slice
	arma::mat mat_grad_precMat  = arma::sum(grad_precMat_cube, 0);  //dim=0 indicates the sum of elements in each column within each slice
	grad_mu				= mat_grad_m.t(); //sum over proteins = dim 0 --> transpose --> matrix(comp, 400)
	grad_precMat		= mat_grad_precMat.t(); //sum over proteins = dim 0 --> transpose --> matrix(comp, 400)

	/*
	 * Check if summation above uses correct dimensions!
	 *
	 *double sum_w =0;
	 *double double sum_m =0;
	 *double double sum_v =0;
		for(int p =0; p < nr_proteins; p++){
			sum_w += grad_weight_contact_mat(p,1);
			sum_m += grad_mu_cube(p,1,0);
			sum_v += grad_precMat_cube(p,1,0);
		}
		std::cout  << "sum_w " << sum_w << std::endl;
		std::cout  << "weight " <<grad_weight_contact(1) << std::endl;
		std::cout  << "sum_m " << sum_m << std::endl;
		std::cout  << "mu " << grad_mu(1,0) << std::endl;
		std::cout  << "sum_v " << sum_v << std::endl;
		std::cout  << "var " <<grad_precMat_cube(1,0) << std::endl;
	 */




}






