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
	input_node u = input_node(&basis_poly, 1, "u");
	
	const_node c1 = const_node("c1");
	const_node c2 = const_node("c2");
	const_node c3 = const_node("c3");

	mult_node a = mult_node(&basis_poly, "a");
	delay_node b = delay_node("b");
	mult_node c = mult_node(&basis_poly, "c");
	add_node d = add_node("d");
	add_node y = add_node("y");
	delay_node y_d = delay_node("y_d");
	mult_node y_d_c = mult_node(&basis_poly, "y_d_c");
	mult_node e = mult_node(&basis_poly, "e");

	// Connect/Init Netlist
	c1.init(0.05);
	c2.init(0.25);
	c3.init(0.6);

	a.init(&u, &c1);
	b.init(&u);
	c.init(&b, &c2);
	d.init(&c, &e);
	y.init(&d, &y_d_c);
	y_d.init(&y);
	y_d_c.init(&y_d, &c3);
	e.init(&a, &y_d);

	// Set bitwidths
	a.set_bitwidth(2);
	c.set_bitwidth(1);
	y_d_c.set_bitwidth(1);
	e.set_bitwidth(1);

	// Generate Basis Polynomials
	basis_poly.generate_polys(2);

	// Set up Simulator
	Simulator sim = Simulator(16);
	sim.add_basis_poly_set(basis_poly);

	sim.add_node(u);
	sim.add_node(c1);
	sim.add_node(c2);
	sim.add_node(c3);
	sim.add_node(a);
	sim.add_node(b);
	sim.add_node(c);
	sim.add_node(d);
	sim.add_node(y);
	sim.add_node(y_d);
	sim.add_node(y_d_c);
	sim.add_node(e);

	// Run Sim
	int tot_sim_steps = 80;
	sim.set_sim_params(MONTE_CARLO, tot_sim_steps);
	sim.run_sim(&u);

	sim.set_sim_params(PCE, tot_sim_steps);
	sim.run_sim(&u);

	// sim.add_plot_node(u);
	sim.add_plot_node(y);
	sim.plot();

	return 0;
}