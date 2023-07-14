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
	input_node eps(&basis_poly, 1, "eps");

	const_node c1(&basis_poly, "c1");
	const_node c2(&basis_poly, "c2");

	mult_node a(&basis_poly, "a");
	mult_node b(&basis_poly, "b");
	mult_node c(&basis_poly, "c");
	mult_node d(&basis_poly, "d");

	mult_node e(&basis_poly, "e");
	mult_node f(&basis_poly, "f");
	sub_node g(&basis_poly, "g");

	// Connect/Init Netlist
	c1.init(10.0);
	c2.init(1.0);
	a.init(&eps, &c1);
	b.init(&eps, &c1);
	c.init(&eps, &c2);
	d.init(&eps, &c2);
	e.init(&a, &b);
	f.init(&c, &d);
	g.init(&e, &f);

	// Generate Basis Polynomials
	basis_poly.generate_polys(3);

	//// Set up Simulator
	Simulator sim = Simulator();
	sim.add_basis_poly_set(basis_poly);

	// Add nodes to simulator
	sim.add_node(eps);
	sim.add_node(c1);
	sim.add_node(c2);
	sim.add_node(a);
	sim.add_node(b);
	sim.add_node(c);
	sim.add_node(d);
	sim.add_node(e);
	sim.add_node(f);
	sim.add_node(g);

	// Set bitwidths
	sim.set_bitwidth(eps, 10);
	sim.set_bitwidth(e, 10);
	sim.set_bitwidth(f, 10);
	sim.set_bitwidth(g, 10);

	// Run Sim
	int tot_sim_steps = 10;
	// sim.set_sim_params(MONTE_CARLO, tot_sim_steps);
	// sim.run_sim(&pll_in);

	sim.set_sim_params(PCE, tot_sim_steps);
	sim.run_sim(&eps);

	sim.add_plot_node(*g.tail);
	sim.plot();


	return 0;
}