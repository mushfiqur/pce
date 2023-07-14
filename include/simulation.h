#ifndef _SIMULATION_H
#define _SIMULATION_H

#include <deque>
#include <algorithm>
#include <random>
#include <vector>

#include "../include/enums.h"
#include "../include/dfg_node.h"
#include "../include/BasisPolys.h"

struct bitwidth_config{
	dfg_node* n;
	int bitwidth;
};

class Simulator{
	public:
	std::vector<bitwidth_config> curr_config;

	public:
	Simulator();
	void add_basis_poly_set(BasisPolySet& bp_set);
	void add_node(dfg_node& n);
	void set_bitwidth(dfg_node& n, int bitwidth);
	void set_sim_params(SimType sim_t, int tot_sim_steps, int mc_samples=100000);
	void run_sim(dfg_node* n);
	void add_plot_node(dfg_node& n);
	void plot();

	private:
	int tot_sim_steps;
	int mc_samples;
	SimType sim_t;
	std::vector<dfg_node*> nodes_arr;
	std::vector<dfg_node*> nodes_to_plot;
	BasisPolySet* bp_set_ptr;

	bool mc_sim_done;
	bool pce_sim_done;

	void set_node_sim_params();
	void calc_bitwidths();
	void propagate_coeffs(dfg_node* n);
};

#endif