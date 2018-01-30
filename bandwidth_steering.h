//#include <vector>
#include <stdafx.h>
//#include <optimization.h>
#include <ap.h>
/**
 * A is just 2*I since the minimization multiplies the
 **/
namespace optimization {

struct optimization_parameters{
	// minimize ||Ax - b||
	// subject to Cx <= C[n+1]
	alglib::real_2d_array A;
	alglib::real_1d_array b;
	alglib::real_2d_array C;
	alglib::integer_1d_array ct; // basically tells the optimizer the i-th constrain has what inequality sign
};

void configure_constraints(int num_groups, 
							std::string& C_str, 
							std::string& d_str, 
							optimization_parameters* opt_params,
							std::vector<std::vector<float>>& traffic_matrix);

void configure_parameters(int num_groups, 
							optimization_parameters* opt_params, 
							std::vector< std::vector<float> >& traffic_matrix);
	


void optimize_allocation(int num_groups, 
							optimization_parameters* opt_params, 
							std::vector< std::vector<float> >& traffic_matrix);

};
