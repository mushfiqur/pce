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
	BasisPolySet basis_poly = BasisPolySet(GAUSSIAN);

	input_node x(&basis_poly, 1, "x");
	input_node y(&basis_poly, 1, "y");
	mult_node z(&basis_poly, "z");
	divide_node out(&basis_poly, "out");
	
	// Connect/Init Netlist
	z.init(&x,&y);
	out.init(&z, &z);

	// Generate Basis Polynomials
	basis_poly.generate_polys(2);

	// Set up Simulator
	Simulator sim = Simulator();
	sim.add_basis_poly_set(basis_poly);

	sim.add_node(x);
	sim.add_node(y);
	sim.add_node(z);
	sim.add_node(out);

	// Set bitwidths

	// Set output node
	sim.set_output_node(out);

	// Run Sim
	int tot_sim_steps = 1;

	sim.set_sim_params(PCE, tot_sim_steps);
	sim.run_sim(&x);
	
	// x.print();
	out.print();

	return 0;
}