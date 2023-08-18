#ifndef _DSP_BLOCK_H
#define _DSP_BLOCK_H

#include <vector>
#include <iostream>
#include <string>
#include <random>
#include "gnuplot-iostream.h"

#include "../include/enums.h"
#include "../include/dfg_node.h"
#include "../include/BasisPolys.h"

class sine_node : public dfg_node {
	private:
	std::vector<const_node*> constants;
	std::vector<dfg_node*> nodes;

	public:
	int a;
	sine_node(BasisPolySet* bp_set, std::string label);
	void init(dfg_node* arg_node);
	void process(int curr_timestamp) override;
	void set_bitwidth(int width) override;
	void set_sim_params(int tot_sim_steps, int mc_samples, int basis_set_size, SimType sim_type) override;
	void print(bool print_last=false);

	void reorder_signal_polys() override;
	void save_signal_polys() override;
	void remove_signal_component() override;

	double get_pwr() override;

	~sine_node();
};


class cosine_node : public dfg_node {
	private:
	std::vector<const_node*> constants;
	std::vector<dfg_node*> nodes;

	public:
	int a;
	cosine_node(BasisPolySet* bp_set, std::string label);
	void init(dfg_node* arg_node);
	void process(int curr_timestamp) override;
	void set_bitwidth(int width) override;
	void set_sim_params(int tot_sim_steps, int mc_samples, int basis_set_size, SimType sim_type) override;
	void print(bool print_last=false);

	void reorder_signal_polys() override;
	void save_signal_polys() override;
	void remove_signal_component() override;

	double get_pwr() override;
	void get_pce_stats(std::vector<double>& mean_arr, std::vector<double>& var_arr) override;


	~cosine_node();
};


#endif