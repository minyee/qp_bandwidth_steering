#include "connectivity_matrix.h"
#include <iostream>
#include <fstream>
#include <cassert>
#include <climits>
//#include "bandwidth_steering.h"

using namespace std;

/**
 * searches through an array with size @size and finds the minimum element and returns the index 
 * of the minimum value, and stores the acutal value of the minimum in solution.
 **/
int64_t search_min(int64_t array[], int64_t size, int64_t& solution) {
	if (size == 1) {
		solution = array[0];
		return 0;
	}

	int64_t index = -1;
	int64_t current_min = INT_MAX;
	for (int64_t i = 0; i < size; i++) {
		if (array[i] == -1)
			continue;
		if (current_min > array[i] && array[i] != -1) {
			index = i;
			current_min = array[i];
		}
	}
	array[index] = -1;
	return index;
}

void target_groups(int64_t num_groups, 
					int64_t src_group, 
					std::vector<int64_t>& incoming_switches_left, 
					std::vector<int64_t>& groups) {

	int64_t target_groups[num_groups-1];
	int64_t index = 0;
	// given a src_group, forms an array of all
	for (int64_t i = 0; i < num_groups; i++) {
		if (src_group == i) 
			continue;
		target_groups[index] = i;
		index++;
	}


	int64_t size = num_groups - 1;
	int64_t array[size];
	for (int64_t i = 0; i < size; i++) {
		array[i] = incoming_switches_left[target_groups[i]];
	}
	for (int64_t i = 0; i < size; i++) {
		int64_t solution;
		int64_t index = search_min(array, size, solution);
		groups.push_back(target_groups[index]);
	}
}

void generate_butterfly(int64_t num_groups,  std::vector<node*>& groups, std::vector<node*>& optical_switches) {
	int64_t num_optical_switches = num_groups - 1;
	for (int64_t i = 0; i < num_groups; i++) {
		groups[i] = new node(i, num_optical_switches, false);
		assert(groups[i]);
	}
	for (int64_t i = 0; i < num_optical_switches; i++) {
		optical_switches[i] = new node(i, num_groups, true);
		assert(optical_switches[i]);
	}
};

void link_butterfly(int64_t num_groups, std::vector<node*>& groups, std::vector<node*>& optical_switches) {
	for (int64_t i = 0; i < num_groups; i++) {
		node* curr_group_node = groups[i];
		assert(curr_group_node);
		assert(!curr_group_node->is_optical_switch());
		for (int64_t j = 0; j < num_groups - 1; j++) {
			node* curr_optical_switch = optical_switches[j];
			assert(optical_switches[j]);
			curr_group_node->set_child(j, curr_optical_switch);
			printf(" The point64_ter for optical_switch: %lld is %p\n", j, curr_optical_switch);
		}
	}
	for (int64_t i = 0; i < num_groups - 1; i++) {
		node* current_optical_switch = optical_switches[i];
		assert(current_optical_switch);
		assert(current_optical_switch->is_optical_switch());
		for (int64_t j = 0; j < num_groups; j++) {
			node* current_group = groups[j];
			assert(current_group);
			printf(" The point64_ter for group: %lld is %p\n", j, current_group);
			current_optical_switch->set_child(j, current_group);
		}
	}
};

/**
 * Butterfly is the canonical Dragonfly case in which every group should be able to get 
 * to every other group in the int64_terconnect using just 1 int64_ter-group link. This function 
 * both generates the group nodes and the optical switch nodes and then subsequently links 
 * them together. 
 * NOTE: This function has to always be followed up using a delete_buttefly function
 * 			call in order to free all the nodes which were malloc-ed earlier using form butterfly
 **/
void form_butterfly(int64_t num_groups, std::vector<node*>& groups, std::vector<node*>& optical_switches) {
	generate_butterfly(num_groups, groups, optical_switches);
	link_butterfly(num_groups, groups,optical_switches);
	for (int64_t i = 0; i < num_groups; i++) {
		node* group = groups[i];
		assert(!group->is_optical_switch());
		std::cout << "current group is: " << group->get_id() << std::endl;
		for (int64_t j = 0; j < num_groups - 1; j++) {
			printf("port: %lld is connected to switch: %lld with point64_ter %p\n",j,optical_switches[j]->get_id(),optical_switches[j]);
			//std::cout << "port: " << std::to_string(j) << " connects to switch: " << std::to_string((group->get_child(j))->get_id()) << std::endl;
		} 
	}

	for (int64_t i = 0; i < num_groups - 1; i++) {
		node* optical_switch = optical_switches[i];
		std::cout << "current optical_switch is: " << optical_switch->get_id() << std::endl;
		assert(optical_switch->is_optical_switch());
		for (int64_t j = 0; j < num_groups; j++) {
			printf("port: %lld is connected to group: %lld with point64_ter %p\n",j,groups[j]->get_id(),groups[j]);
		}
	}
	//for (int64_t i)
};

void delete_butterfly(int64_t num_groups, std::vector<node*>& groups, std::vector<node*>& optical_switches) {
	for (int64_t i = 0; i < groups.size(); i++) {
		delete groups[i];
	}	
	for (int64_t i = 0; i < optical_switches.size(); i++) {
		delete optical_switches[i];
	}
};

/**
 * The algorithm that calls the DFS function to search the butterfly graph 
 * that is the Flexfly topology to figure out what the configuration should 
 * be 
 **/
void canonical_dragonfly_config_greedy(int64_t num_groups, std::vector<node*>& groups, 
										std::vector<node*>& optical_switches,  
										std::vector<std::vector<int64_t> >& optical_inout_connections) {
	// for all switches, find the optical switch that gets them to every other switch
	std::vector<int64_t> incoming_switches_left (num_groups);
	std::fill(incoming_switches_left.begin(), incoming_switches_left.end(), num_groups - 1);
	for (int64_t i = 0; i < num_groups; i++) {
		node* curr_group = groups[i];
		assert(curr_group);
		int64_t curr_group_id = curr_group->get_id();
		std::vector<int64_t> in_order_group;
		
		// prioritize the 
		target_groups(num_groups, i, incoming_switches_left, in_order_group);
		for (int64_t j = 0; j < num_groups - 1; j++) {	
			int64_t target_group = in_order_group[j];
			// NOTE: j is the target group's id
			node* optical_switch = dfs(curr_group, 3, target_group);
			assert(curr_group->get_id() != target_group);
			if (optical_switch == nullptr) {
				std::cout << "the failed switch: " << std::to_string(curr_group->get_id()) << " and target group: " << std::to_string(target_group) << std::endl;
			} else {
				incoming_switches_left[target_group]--; // NOTE: keep this in mind
				assert(optical_switch->is_optical_switch());
				std::vector<int64_t>& inout_connection_for_optical_switch = optical_inout_connections[optical_switch->get_id()];
				if (inout_connection_for_optical_switch.size() != num_groups) {
					inout_connection_for_optical_switch.resize(num_groups);
				}
				inout_connection_for_optical_switch[i] = target_group; 
			}
		}
	}
};

/**
 * short for depth-first search
 */ 
node* dfs(node* curr_node, int64_t depth, int64_t target_group_id) {
	if (curr_node == nullptr) 
		return nullptr;
	
	node* found_node = nullptr;
	if (depth == 1) {
		assert(!curr_node->is_optical_switch());
		found_node = (target_group_id == curr_node->get_id()) ? curr_node : nullptr;
	} else {
		//cout << to_string(curr_node->get_num_child()) << endl;
		int64_t iter = curr_node->get_num_child();
		for (int64_t i = 0; i < iter; i++) {
			found_node = dfs(curr_node->get_child(i), depth - 1, target_group_id);
			if (found_node == nullptr)
				continue;
			
				// return curr node so that the group will know which 
				// optical switch can be connected to the destination group
			curr_node->set_child(i, nullptr);
			if (curr_node->is_optical_switch()) {
				assert(depth == 2);
				found_node = curr_node;
			}
			break;
		}
	}
	return found_node;
};


void check_switch_config(std::vector< std::vector <int64_t> >& optical_inout_connections) {
	for (int64_t i = 0; i < optical_inout_connections.size(); i++) {
		std::vector<int64_t>& current_config = optical_inout_connections[i];
		std::cout << "Switch: " << std::to_string(i) << std::endl;
		for (int64_t j = 0; j < current_config.size(); j++) {
			std::cout << "inport: " << std::to_string(j) << " -> outport: " << std::to_string(current_config[j]) << std::endl; 
		}
		std::cout << std::endl;
	}	
}

/**
 * Combines all the above the algorithms above and int64_to one. It begins by generating the physical
 * Butterfly topology connection where middle layer switches are all optical-switches. Next it will
 * run a DFS algorithm to find the canonical configuration (dfly with 1 int64_ter-group link to every group)
 * and finally frees all the nodes that were allocated when calling form_butterfly() function.
 *
 * NOTE: This is the function that should really be called by external codes as it abstracts away a lot of 
 * 			the int64_ternal implementation of finding the optical switch inport-outport configuration needed to
 * 			achieve the canonical dragonfly topology.
 **/
void configure_optical_switches_canonical(int64_t num_groups, std::vector< std::vector<int64_t> >& optical_inout_connections) {
	std::vector<node*> groups;
	std::vector<node*> optical_switches;
	groups.resize(num_groups);
	groups.reserve(num_groups);
	optical_switches.resize(num_groups - 1);
	optical_switches.reserve(num_groups - 1);
	// first form (generate and link) butterfly physical connection
	form_butterfly(num_groups, groups, optical_switches);
	optical_inout_connections.resize(num_groups - 1);
	optical_inout_connections.reserve(num_groups - 1);
	// traverse the butterfly graph to find the canonical dfly optical configuration
	canonical_dragonfly_config_greedy(num_groups, groups, optical_switches, optical_inout_connections);
	// finally delete all the groups and optical switch nodes;
	delete_butterfly(num_groups, groups, optical_switches);
}

inline int64_t wrapped_increment(int64_t index, int64_t max_num) {
	return (index + 1) % max_num;
}

void configure_simpler_model(int64_t num_groups, std::vector< std::vector<int64_t> >& optical_inout_connections) {
	int64_t num_optical_switches = num_groups - 1;
	if (optical_inout_connections.size() == 0) {
		//optical_inout_connections.reserve(num_optical_switches);
		optical_inout_connections.resize(num_optical_switches);
		optical_inout_connections.reserve(num_optical_switches);
	}
	for (int64_t i = 0; i < num_groups; i++) {
		// i equals the current group
		int64_t tmp = i; // stores the current index
		for (int64_t j = 0; j < num_groups - 1; j++) {
			std::vector<int64_t>& current_switch_config = optical_inout_connections[j];
			if (current_switch_config.size() == 0) {
				//current_switch_config.reserve(num_groups);
				current_switch_config.resize(num_groups);
				current_switch_config.reserve(num_groups);
			}
			tmp = wrapped_increment(tmp, num_groups);
			current_switch_config[i] = tmp;
		}
		
	}
}

void form_traffic_matrix(int64_t num_groups, std::vector< std::vector<int64_t> >& traffic_matrix) {
	//traffic_matrix.resize(num_groups);
	for (int64_t i = 0; i < num_groups; i++) {
		std::vector<int64_t> tmp;
		//tmp.resize(num_groups);
		for (int64_t j = 0; j < num_groups; j++) {
			if (i == j) {
				tmp.push_back(0);
			} else {
				tmp.push_back(1);
				//tmp.push_back(0);
			}
		}
		traffic_matrix.push_back(tmp);
	}
};

void read_ETM_from_file(std::string filename, std::vector< std::vector<double> >& traffic_matrix, int64_t& num_groups) {
	ifstream file_stream (filename);
	// file stream has to be opened
	if (!file_stream.is_open()) {
		std::cout << "file cannot be opened" << std::endl;
	}
	int64_t matrix_size;
	file_stream >> matrix_size; 
	num_groups = matrix_size;
	int64_t entry;
	int64_t index = 0;
	int64_t row = 0;
	int64_t col = 0;
	std::cout << "The matrix size is: " << std::to_string(matrix_size) << std::endl;
	for (int64_t i = 0; i < matrix_size; i++) {
		std::vector<double> row_entries;
		for (int64_t j = 0; j < matrix_size; j++) {
			row_entries.push_back(0);
		}
		traffic_matrix.push_back(row_entries);
	}
	//file_stream >> matrix_size;
	//std::cout << "testing second entry " << std::to_string(matrix_size) << std::endl;
	while (!file_stream.eof()) {
		std::cout << "Current iteration: " << std::to_string(index) << std::endl;
		file_stream >> entry;
		std::cout << "Current entry: " << std::to_string(entry) << std::endl;
		row = index / matrix_size;
		col = index % matrix_size;
		std::vector<double>& row_entries = traffic_matrix[row];
		row_entries[col] = entry;
		index++;
	}

	cout << "index is: " << to_string(index) << endl;
	assert(index == (matrix_size * matrix_size));
	file_stream.close();
};

void massage_matrix(std::vector<std::vector<double>>& traffic_matrix) {
	for (int64_t i = 0; i < traffic_matrix.size(); i++) {
		double sum = 0;
		for (int64_t j = 0; j < traffic_matrix[i].size(); j++) {
			sum += traffic_matrix[i][j];
		} 
		for (int64_t j = 0; j < traffic_matrix[i].size(); j++) {
			traffic_matrix[i][j] = (traffic_matrix[i][j] / sum) * ((traffic_matrix.size() - 1));
		}
	}
	cout << "Print64_ting matrix" << endl;
	for (auto row : traffic_matrix) {
		cout << endl;
		for (auto entry : row) {
			cout << to_string(entry) << " ";
		}
	}
};
//void round_down(std::vector<std::vector<>>)
/**
 * Tests the correctness of the DFS algorithm that is used to search for a connection pattern 
 * 
 **/
int main(int argc, char** argv) {
	//if (argc != 2) {
	//	cout<< "usage: "<< argv[0] <<" <filename>\n";
	//}
	int64_t massage;
	if (argc > 2) {
		massage = stoi(argv[2]);
	}
	int64_t num_groups;
	/*
	vector< vector<int64_t> > optical_inout_connections;
	configure_optical_switches_canonical(num_groups, optical_inout_connections);
	check_switch_config(optical_inout_connections);
	vector< vector<int64_t> > optical_inout_connections1;
	configure_simpler_model(num_groups, optical_inout_connections1);
	check_switch_config(optical_inout_connections1);
	std::cout << "Exited perfectly, good job!" << std::endl;
	*/
	std::vector<std::vector<double>> traffic_matrix;
	read_ETM_from_file(argv[1], traffic_matrix, num_groups);
	if (massage > 0) {
		massage_matrix(traffic_matrix);
	}
	//cout << "hello 2 " << endl;
	
	optimization::optimization_parameters opt_params;
	std::string C_str;
	std::string d_str;
	configure_parameters(num_groups, &opt_params, traffic_matrix);
	configure_constraints(num_groups, C_str, d_str, &opt_params,traffic_matrix);
	std::cout << C_str << std::endl;
	std::cout << d_str << std::endl;
	//std::cout << opt_params.b.tostring(2) << std::endl;
	optimization::optimize_allocation(num_groups, &opt_params, traffic_matrix);

	return 0;
}
