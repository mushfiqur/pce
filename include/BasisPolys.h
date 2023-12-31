#ifndef _BASIS_POLYS_H
#define _BASIS_POLYS_H

#include <iostream>
#include <vector>
#include <string>
#include <math.h>

#include "../include/enums.h"
#include "../include/polynomial.h"

class BasisPolySet {
	public:
	std::vector<var*> var_arr;
	std::vector<var*> noise_vars;
	std::vector<polynomial*> basis_polys;

	std::vector<double> poly_sqr_expt;
	std::vector<double> poly_expt;
	std::vector<std::vector<std::vector<double>>> expt_table;

	int num_vars;
	int order;
	int set_size;

	BasisPolySet(RandVarDist dist_type);
	~BasisPolySet();

	void add_variable(var* v);
	void generate_polys(int order);
	void regenerate_polys(int order);
	
	int get_var_idx(int id);
	int get_var_idx(int id1, int id2);
	int get_new_var_id();
	
	void print();
	void print_polys();
	void print_exp_table();

	private:
	int max_var_idx;

	void gen_exp_table();
	void gen_poly_expt_sqr_table();
	void gen_poly_expt_table();

};

#endif