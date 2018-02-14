#include "data_structures.h"
#include <vector>
#include "bandwidth_steering.h"
void generate_butterfly(int64_t num_groups,  std::vector<node*>& groups, std::vector<node*>& optical_switches);

void link_butterfly(int64_t num_groups, std::vector<node*>& groups, std::vector<node*>& optical_switches);

void form_butterfly(int64_t num_groups, std::vector<node*>& groups, std::vector<node*>& optical_switches);

void delete_butterfly(int64_t num_groups, std::vector<node*>& groups, std::vector<node*>& optical_switches);

void canonical_dragonfly_config_greedy(int64_t num_groups, std::vector<node*>& groups, 
										std::vector<node*>& optical_switches,  
										std::vector<std::vector<int64_t> >& optical_inout_connections);

void configure_optical_switches_canonical(int64_t num_groups, std::vector< std::vector<int64_t> >& optical_inout_connections);

node* dfs(node* curr_node, int64_t depth, int64_t target_group_id);

void massage_matrix(std::vector<std::vector<double>>& traffic_matrix);
