#ifndef _BASIS_POLYS_H
#define _BASIS_POLYS_H

#include <iostream>
#include <vector>
#include <string>
#include <math.h>

#include "../include/polynomial.h"


// #define USE_HERMITE_POLYS
#define USE_LEGENDRE_POLYS

class BasisPolySet {
	public:
	std::vector<var*> var_arr;
	std::vector<polynomial*> basis_polys;

	std::vector<double> poly_sqr_expt;
	std::vector<double> poly_expt;
	std::vector<std::vector<std::vector<double>>> expt_table;

	int num_vars;
	int order;
	int set_size;

	BasisPolySet(int num_vars, int order);
	~BasisPolySet();
	
	int get_var_idx(int id);
	
	void print();
	void print_polys();
	void print_exp_table();

	private:
	void gen_exp_table();
	void gen_poly_expt_sqr_table();
	void gen_poly_expt_table();

};

#endif