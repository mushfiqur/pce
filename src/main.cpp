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
	// Define Netlist
	BasisPolySet basis_poly = BasisPolySet(UNIFORM);

	input_node x(&basis_poly, 1, "x");
	mult_node x_c(&basis_poly, "x_c");
	const_node c1(&basis_poly, "c1");
	const_node c2(&basis_poly, "c2");
	
	delay_node out_d(&basis_poly, "y_d");
	mult_node out_dc(&basis_poly, "y_dc");
	add_node out(&basis_poly, "out");

	// Connect/Init Netlist
	c1.init(0.5);
	c2.init(0.5);

	x_c.init(&c2, &x);
	out_d.init(&out);
	out_dc.init(&c1, &out_d);
	out.init(&out_dc, &x_c);

	// Generate Basis Polynomials
	basis_poly.generate_polys(3);

	// Set up Simulator
	Simulator sim = Simulator();
	sim.add_basis_poly_set(basis_poly);

	sim.add_node(x);
	sim.add_node(x_c);
	sim.add_node(c1);
	sim.add_node(c2);
	sim.add_node(out_d);
	sim.add_node(out_dc);
	sim.add_node(out);

	// Set bitwidths
	sim.set_bitwidth(x, 32);
	sim.set_bitwidth(x_c, 32);
	sim.set_bitwidth(out_dc, 32);

	sim.set_output_node(out);

	// Run Sim
	int tot_sim_steps = 10;

	sim.set_sim_params(PCE, tot_sim_steps);
	// sim.run_sim(&b);
	sim.run_sim_anneal(&x, 1.0/3.0, 20.0, 100000);
	
	sim.print();

	// sim.add_plot_node(out);
	// sim.plot();

	return 0;
}