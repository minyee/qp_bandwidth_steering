#include "bandwidth_steering.h"
#include <optimization.h>
//#include <stdafx.h>
//#include <ap.h>

namespace optimization {

std::vector<double> vector_sum(std::vector<double>& a, std::vector<double>& b, int64_t size){
	std::vector<double> sol;
	for (int64_t i = 0; i < size; i++) {
		sol.push_back(a[i] + b[i]);
	}
	return sol;
};

// calculates a transposed multiplied by b, a.k.a transpose(a) * b.
int64_t inner_product(std::vector<int64_t>& a, std::vector<int64_t>& b, int64_t size) {
	int64_t sol = 0;
	for (int64_t i = 0; i < size; i++) {
		sol += a[i] * b[i];
	}
	return sol;
};

std::vector<int64_t> matrix_vector_multiply(std::vector< std::vector<int64_t> >& matrix, std::vector<int64_t>& a, int64_t num_rows, int64_t num_cols) {
	std::vector<int64_t> sol;
	for (int64_t i = 0; i < num_rows; i++) {
		sol.push_back(inner_product(matrix[i],a,num_cols));
	}
	return sol;
};

inline bool row_inrange(int64_t index, int64_t group, int64_t num_groups) {
	return index >= group * num_groups && index < (group + 1) * num_groups ? true : false;
};	

void configure_constraint64_ts_individual(int64_t num_groups, 
								optimization_parameters* opt_params, 
								std::string& C_str, 
								std::string& d_str,
								std::vector< std::vector<double> >& traffic_matrix) {
	int64_t n = num_groups * num_groups;
	for (int64_t i = 0 ; i < n; i++) {
		C_str += "[";
		for (int64_t j = 0; j < n; j++) {
			if (i == j) {
				C_str += "1,";
			} else {
				C_str += "0,";
			} 
		}
		C_str += std::to_string((num_groups - 1) - traffic_matrix[i/num_groups][i%num_groups]);
		C_str += "],";
		C_str += "[";
		d_str += "-1,";
		for (int64_t j = 0; j < n; j++) {
			if (i == j) {
				C_str += "1,";
			} else {
				C_str += "0,";
			} 
		}
		C_str += std::to_string(-traffic_matrix[i/num_groups][i%num_groups]);
		C_str += "],";
		d_str += "1,";
	}
}

/**
 * num_groups is the number of groups in the topology
 * This happens to be the same as the
 * NOTE: THIS HAS TO BE CALLED BEFORE configure_constraint64_ts_cols
 **/
void configure_constraint64_ts_rows(int64_t num_groups, 
								optimization_parameters* opt_params, 
								std::string& C_str, 
								std::string& d_str) {
	if (!opt_params)
		return;
	int64_t n = num_groups * num_groups;
	for (int64_t i = 0; i < num_groups; i++) {
		C_str += "[";
		for (int64_t index = 0; index < n; index++) {
			if (row_inrange(index, i, num_groups)) {
				C_str += "1";
			} else {
				C_str += "0";
			} 
			C_str += ",";
		}
		//uses the invariant that for each traffic matrix's row will ultimately sum to num_groups - 1
		C_str += "0";
		//C_str += std::to_string(2*(num_groups - 1));
		//C_str += std::to_string(num_groups - 1);
		d_str += "-1,";
		C_str += "],\n";
	}

	for (int64_t i = 0; i < num_groups; i++) {
		C_str += "[";
		int64_t row_start = i * n;
		for (int64_t index = 0; index < n; index++) {
			if (row_inrange(index, i, num_groups)) {
				C_str += "1";
			} else {
				C_str += "0";
			} 
			C_str += ",";
		}
		//uses the invariant that for each traffic matrix's row will ultimately sum to num_groups - 1
		//C_str += "-";
		C_str += std::to_string(-(num_groups - 1));
		//C_str += std::to_string(2*(num_groups - 1));
		//C_str += std::to_string(num_groups - 1);
		d_str += "1,";
		C_str += "],\n";
	}
};

void configure_constraint64_ts_cols(int64_t num_groups, 
								optimization_parameters* opt_params, 
								std::string& C_str, 
								std::string& d_str,
								std::vector< std::vector<double> >& traffic_matrix) {
	if (!opt_params)
		return;
	int64_t n = num_groups * num_groups;
	
	for (int64_t col = 0; col < num_groups; col++) {
		int64_t cnt = 0;
		C_str += "[";
		for (int64_t i = 0; i < n; i++) {
			if ((cnt * num_groups + col) == i) {
				C_str += "1,";
				cnt++;
			} else {
				C_str += "0,";
			}
		}
		int64_t bound = num_groups - 1;
		for (int64_t row = 0; row < num_groups; row++) {
			bound -= traffic_matrix[row][col];
		}
		C_str += std::to_string(bound);
		cnt = 0;
		if (col < num_groups - 1) {
			C_str += "],\n";
			
		} else {
			C_str += "],\n";
			//d_str += "-1";
		}
		d_str += "-1,";
		
	}	
	for (int64_t col = 0; col < num_groups; col++) {
		int64_t cnt = 0;
		C_str += "[";
		for (int64_t i = 0; i < n; i++) {
			if ((cnt * num_groups + col) == i) {
				C_str += "1,";
				cnt++;
			} else {
				C_str += "0,";
			}
		}
		int64_t bound = 0;
		for (int64_t row = 0; row < num_groups; row++) {
			bound -= traffic_matrix[row][col];
		}
		C_str += std::to_string(bound);
		//std::cout << std::to_string(bound) << std::endl;
		cnt = 0;
		if (col < num_groups - 1) {
			C_str += "],\n";
			d_str += "1,";
		}
		else {
			C_str += "]";
			d_str += "1";
		}
		
	}
};

void configure_constraints(int64_t num_groups, 
							std::string& C_str, 
							std::string& d_str, 
							optimization_parameters* opt_params, 
							std::vector<std::vector<double>>& traffic_matrix) {
	C_str += "[";
	d_str += "[";
	configure_constraint64_ts_individual(num_groups, opt_params, C_str, d_str, traffic_matrix);
	configure_constraint64_ts_rows(num_groups, opt_params, C_str, d_str);
	configure_constraint64_ts_cols(num_groups, opt_params, C_str, d_str, traffic_matrix);
	int64_t num_constraint64_ts = 2 * num_groups;
	/*
	for (int64_t i = 0; i < num_constraint64_ts; i++) {
		if (i < num_constraint64_ts - 1) 
			d_str += "-1,";
		else 
			d_str += "-1";

	}
	*/
	C_str += "]";
	d_str += "]";
	//std::cout << C_str << std::endl;
	//std::cout << d_str << std::endl;
};

void configure_parameters(int64_t num_groups, optimization_parameters* opt_params, std::vector< std::vector<double> >& traffic_matrix) {
	if (!opt_params) 
		return;
	//opt_params->C = "";
	int64_t n = num_groups * num_groups;
	std::string A_str = "[";
	std::string b_str = "[";
	for (int64_t i = 0; i < n; i++) {
		//std::vector<int64_t>& traffic_row = traffic_matrix[i];
		int64_t row = i / num_groups;
		int64_t col = i % num_groups;
		A_str += "[";
		b_str += "0";
		//b_str += std::to_string(-(traffic_matrix[row])[col]);
		if (i < n - 1)
			b_str += ",";
		for (int64_t j = 0; j < n; j++) {
			if (i==j)
				A_str += "2";
			else
				A_str += "0";
			
			if (j < n - 1)
				A_str += ",";
		}
		A_str += "]";
		if (i < n - 1) {
			A_str += ",";
		}
	}
	A_str += "]";
	b_str += "]";
	opt_params->A = alglib::real_2d_array(A_str.c_str());
	std::string C_str;
	std::string d_str;
	configure_constraints(num_groups, C_str, d_str, opt_params, traffic_matrix);
	opt_params->C = alglib::real_2d_array(C_str.c_str());
	opt_params->b = alglib::real_1d_array(b_str.c_str());
	opt_params->ct = alglib::integer_1d_array(d_str.c_str());
};

void optimize_allocation(int64_t num_groups, optimization_parameters* opt_params, std::vector< std::vector<double> >& traffic_matrix) {
	alglib::minqpstate state;
    alglib::minqpreport rep;
    alglib::real_1d_array x;
    // create solver, set quadratic/linear terms
    alglib::minqpcreate(num_groups*num_groups, state);
    alglib::minqpsetquadraticterm(state, opt_params->A);
    alglib::minqpsetlinearterm(state, opt_params->b);
    alglib::minqpsetlc(state, opt_params->C, opt_params->ct);

    // Set scale of the parameters.
    // It is strongly recommended that you set scale of your variables.
    // Knowing their scales is essential for evaluation of stopping criteria
    // and for preconditioning of the algorithm steps.
    // You can find more information on scaling at http://www.alglib.net/optimization/scaling.php
    //alglib::minqpsetscale(state, s);

    //
    // Solve problem with BLEIC-based QP solver.
    //
    // This solver is int64_tended for problems with moderate (up to 50) number
    // of general linear constraint64_ts and unlimited number of box constraint64_ts.
    //
    // Default stopping criteria are used.
    //
    alglib::minqpsetalgobleic(state, 0.0, 0.0, 0.0, 0);
    alglib::minqpoptimize(state);
    alglib::minqpresults(state, x, rep);
    std::vector<double> traffic_vector;
    const double* solutions = x.getcontent();
    std::vector<double> soln;
    
    for (int64_t i = 0; i < (num_groups * num_groups); i++) {
    	soln.push_back(solutions[i]);
    	//print64_tf("%f,", soln[i]);
    	int64_t row = i / num_groups;
    	int64_t col = i % num_groups;
    	traffic_vector.push_back(traffic_matrix[row][col]);
    }
    
    std::vector<double> actual_vector = optimization::vector_sum(soln,traffic_vector,(num_groups*num_groups));
    for (int64_t i = 0; i < (num_groups*num_groups); i++) {
    	if (i % num_groups == 0) {
    		printf("]\n[%f,", actual_vector[i]);	
    	//	print64_tf("]\n[%f,", soln[i]);	
    	} else {
    		printf("%f,", actual_vector[i]);	
    		//print64_tf("%f,", soln[i]);	
    	}
    	
    }
    std::cout << "]" << std::endl;
    //std::cout << "The results are: " << x.tostring(2) << std::endl; 
};

}