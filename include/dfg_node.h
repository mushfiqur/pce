#ifndef _DFG_NODE_H
#define _DFG_NODE_H

#include <vector>
#include <iostream>
#include <string>
#include <random>
#include <unsupported/Eigen/IterativeSolvers>
#include "gnuplot-iostream.h"

#include "../include/enums.h"
#include "../include/BasisPolys.h"

typedef enum
{
	ADD,
	SUB,
	MULT,
	DIVIDE,
	INPUT_SIGNAL,
	INPUT_NOISE,
	DELAY,
	CONST
} NodeType;

struct signal_polys {
	double coeff;
	polynomial* poly;
};

class dfg_node{
	public:
	NodeType t;
	SimType sim_type;

	std::string label;
	int last_exec_time;
	int node_id;

	int bitwidth;

	std::vector<std::vector<double>> mc_samples;
	std::vector<std::vector<double>> pce_coeffs;
	std::vector<std::vector<double>> signal_coeffs;
	std::vector<std::vector<signal_polys>> signal;

	std::vector<double> pwr_arr;

	std::vector<dfg_node*> next_nodes;
	std::vector<dfg_node*> prev_nodes;

	dfg_node(NodeType t, BasisPolySet* bp_set, std::string label, int node_id=-1);
	~dfg_node();

	void add_next_node(dfg_node* n);
	void print();
	void print_signal_coeffs();
	void remove_signal_component();
	bool node_args_ready(int curr_timestamp);
	void print_pwr(BasisPolySet& bp_set);
	double get_pwr();
	void get_mc_stats(std::vector<double>& mean_arr, std::vector<double>& var_arr);
	void get_pce_stats(std::vector<double>& mean_arr, std::vector<double>& var_arr, BasisPolySet& bp_set);
	void save_signal_polys();
	void reorder_signal_polys();

	virtual void set_sim_params(int tot_sim_steps, int mc_samples, int basis_set_size, SimType sim_type);
	
	virtual void process(int curr_timestamp);
	virtual void init();

	virtual void set_bitwidth(int width=-1);

	dfg_node* lhs;
	dfg_node* rhs;

	dfg_node* head;
	dfg_node* tail;

	protected:
	BasisPolySet* bp_set_ptr;
	int basis_set_size;
};

class add_node : public dfg_node{
	private:
	void process_pce_sim(int curr_timestamp);
	void process_mc_sim(int curr_timestamp);

	public:
	add_node(BasisPolySet* bp_set, std::string label) : dfg_node(ADD, bp_set, label){};
	void init(dfg_node* lhs, dfg_node* rhs);
	void process(int curr_timestamp) override;
	void set_bitwidth(int width) override;
};

class sub_node : public dfg_node{
	private:
	void process_pce_sim(int curr_timestamp);
	void process_mc_sim(int curr_timestamp);

	public:
	sub_node(BasisPolySet* bp_set, std::string label) : dfg_node(SUB, bp_set, label){};
	void init(dfg_node* lhs, dfg_node* rhs);
	void process(int curr_timestamp) override;
	void set_bitwidth(int width) override;
};

class mult_node : public dfg_node{
	private:
	// std::vector<std::vector<std::vector<double>>>* ref_exp_table;
	// std::vector<double>*  ref_exp_sqr_table;
	
	private:
	void process_pce_sim(int curr_timestamp);
	void process_mc_sim(int curr_timestamp);
	
	public:
	mult_node(BasisPolySet* bp_set, std::string label);
	// ~mult_node();
	void init(dfg_node* lhs, dfg_node* rhs);
	void process(int curr_timestamp) override;
	void set_bitwidth(int width) override;
	void set_sim_params(int tot_sim_steps, int mc_samples, int basis_set_size, SimType sim_type) override;

};

class divide_node : public dfg_node{
	public:
	divide_node(BasisPolySet* bp_set, std::string label);
	void init(dfg_node* lhs, dfg_node* rhs);
	void process(int curr_timestamp) override;
	void set_bitwidth(int width) override;

	private:
	void process_pce_sim(int curr_timestamp);
	void process_mc_sim(int curr_timestamp);

};

class input_node : public dfg_node{
	private:
	var* v;
	double unfm_dist_param_a;
	double unfm_dist_param_b;
	BasisPolySet* bp_set_ptr;

	std::mt19937 mt;
	std::uniform_real_distribution<double> dist;

	public:
	input_node(BasisPolySet* bp_set, int node_id, std::string label);
	void set_range(double a, double b);
	void set_sim_params(int tot_sim_steps, int mc_samples, int basis_set_size, SimType sim_type) override;
	void process(int curr_timestamp) override;
	void init();
	void set_bitwidth(int width) override;
};

class noise_node : public dfg_node{
	private:
	var* v;
	BasisPolySet* bp_set_ptr;

	std::mt19937 mt;
	std::uniform_real_distribution<double> dist;

	public:
	noise_node(BasisPolySet* bp_set, int node_id, std::string label);
	void set_sim_params(int tot_sim_steps, int mc_samples, int basis_set_size, SimType sim_type) override;
	void process(int curr_timestamp) override;
	void init();
	void set_bitwidth(int width) override;
};

class delay_node : public dfg_node{
	private:
	void process_pce_sim(int curr_timestamp);
	void process_mc_sim(int curr_timestamp);

	public:
	delay_node(BasisPolySet* bp_set, std::string label) : dfg_node(DELAY, bp_set, label){};
	void init(dfg_node* prev_node);
	void process(int curr_timestamp) override;
	void set_bitwidth(int width) override;
};

class const_node : public dfg_node{
	public:
	double val;
	const_node(BasisPolySet* bp_set, std::string label) : dfg_node(CONST, bp_set, label){};
	void process(int curr_timestamp) override;
	void init(double val);
	void set_sim_params(int tot_sim_steps, int mc_samples, int basis_set_size, SimType sim_type) override;
};

#endif