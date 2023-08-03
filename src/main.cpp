#include <iostream>
#include <deque>
#include <algorithm>
#include <random>

#include "../include/enums.h"
#include "../include/simulation.h"
#include "../include/utils.h"
#include "../include/BasisPolys.h"
#include "../include/dfg_node.h"

// valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose --log-file=valgrind_log ./main
// valgrind -q --vgdb-error=0 --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose ./main
// valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose ./main > log 2>&1
int main(){
	BasisPolySet basis_poly = BasisPolySet(UNIFORM);
	
	// Define Netlist
	input_node x(&basis_poly, 1, "eps");
	const_node c1(&basis_poly, "c1");

	mult_node x_2(&basis_poly, "x^2");
	mult_node x_3(&basis_poly, "x^3");
	mult_node x_frac(&basis_poly, "x_{frac}");

	sub_node sin_x(&basis_poly, "sin_x");

	divide_node sinc_x(&basis_poly, "sinc_x");

	// Connect/Init Netlist
	c1.init(1.0/6.0);
	x_2.init(&x, &x);
	x_3.init(&x_2, &x);
	x_frac.init(&x_3, &c1);
	sin_x.init(&x, &x_frac);
	sinc_x.init(&sin_x, &x);

	// Generate Basis Polynomials
	basis_poly.generate_polys(2);

	//// Set up Simulator
	Simulator sim = Simulator();
	sim.add_basis_poly_set(basis_poly);

	// Add nodes to simulator
	sim.add_node(x);
	sim.add_node(c1);
	sim.add_node(x_2);
	sim.add_node(x_3);
	sim.add_node(x_frac);
	sim.add_node(sin_x);
	sim.add_node(sinc_x);

	// // Set bitwidths
	// sim.set_bitwidth(x, 4);
	// sim.set_bitwidth(x_2, 4);
	// sim.set_bitwidth(x_3, 4);
	// sim.set_bitwidth(x_frac, 4);

	sim.set_output_node(sinc_x);

	// Run Sim
	int tot_sim_steps = 1;
	// sim.set_sim_params(MONTE_CARLO, tot_sim_steps);
	// sim.run_sim(&pll_in);

	sim.set_sim_params(PCE, tot_sim_steps);
	sim.run_sim(&x);

	sin_x.print();
	x.print();

	sinc_x.print();

	return 0;
}