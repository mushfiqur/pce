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
	double curr_sol_noise_pwr;
	std::vector<bitwidth_config> curr_solution;

	public:
	Simulator();
	void add_basis_poly_set(BasisPolySet& bp_set);
	void add_node(dfg_node& n);
	void set_bitwidth(dfg_node& n, int bitwidth);
	void set_output_node(dfg_node& n);
	void set_sim_params(SimType sim_t, int tot_sim_steps, int mc_samples=100000);
	void run_sim(dfg_node* n);
	void run_sim_anneal(dfg_node* n, double sig_pwr, double tgt_snr, int tot_iters);
	void add_plot_node(dfg_node& n);
	void plot();
	void print();

	private:
	dfg_node* head;
	dfg_node* output_node;

	std::mt19937 mt;
	std::uniform_real_distribution<double> real_dist;
	std::uniform_int_distribution<int> int_dist;

	int tot_sim_steps;
	int mc_samples;
	SimType sim_t;
	std::vector<dfg_node*> nodes_arr;
	std::vector<dfg_node*> nodes_to_plot;
	BasisPolySet* bp_set_ptr;

	bool mc_sim_done;
	bool pce_sim_done;

	std::vector<bitwidth_config> get_neighbour();
	void initialize(dfg_node* n);
	double try_solution(std::vector<bitwidth_config>& proposed_sol);
	void set_node_sim_params();
	void calc_bitwidths();
	void propagate_coeffs();
};

#endif