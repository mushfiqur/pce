#include <iostream>
#include <deque>
#include <algorithm>
#include <random>

#include "../include/enums.h"
#include "../include/simulation.h"
#include "../include/utils.h"
#include "../include/BasisPolys.h"
#include "../include/dfg_node.h"
#include "../include/dsp_block.h"

// valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose --log-file=valgrind_log ./main
// valgrind -q --vgdb-error=0 --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose ./main
// valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose ./main > log 2>&1
#include <iostream>
#include <deque>
#include <algorithm>
#include <random>

#include "../include/enums.h"
#include "../include/simulation.h"
#include "../include/utils.h"
#include "../include/BasisPolys.h"
#include "../include/dfg_node.h"
#include "../include/dsp_block.h"

// valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose --log-file=valgrind_log ./main
// valgrind -q --vgdb-error=0 --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose ./main
// valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose ./main > log 2>&1
int main(){	
	// Define Netlist
	BasisPolySet basis_poly = BasisPolySet(UNIFORM);
	int tot_sim_steps = 1;

	input_node x(&basis_poly, "x");
	fir_node filt(&basis_poly, "filt");
	// const_node c1(&basis_poly, "c1");
	// mult_node dummy(&basis_poly, "dummy");

	// Connect/Init Netlist
	x.add_dist(UNIFORM, tot_sim_steps);
	x.set_range(-1, 1);

	filt.init(&x, "/home/mushf/pce/filters/rf_filt.txt", 1000);
	filt.add_dist(1);

	// c1.init(1.0);
	// dummy.init(&filt, &c1);

	// Generate Basis Polynomials
	basis_poly.generate_polys(2);
	// basis_poly.print();

	// Set up Simulator
	Simulator sim = Simulator();
	sim.add_basis_poly_set(basis_poly);

	sim.add_node(x);
	sim.add_node(filt);
	// sim.add_node(dummy);
	// sim.add_node(c1);


	// // Set bitwidths
	sim.set_bitwidth(x,    32);
	sim.set_bitwidth(filt, 32);
	// sim.set_bitwidth(out_dc, 4);

	sim.set_output_node(filt);

	// Run Sim
	sim.set_sim_params(PCE, tot_sim_steps);
	// sim.run_sim(&x);

	// filt.print();

	sim.run_sim_anneal(&x, 100.0, 10000);

	// out.print_pwr();	
	sim.print();
	// dummy.print_pwr();


	return 0;
}