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
	BasisPolySet basis_poly = BasisPolySet(UNIFORM);
	int tot_sim_steps = 1;
	
	// Define Netlist
	input_node fir_in(&basis_poly, "fir_{in}");

	fir_node filter(&basis_poly, "filter");

	// Connect/Init Netlist
	std::vector<double> coeffs = {0.7640, 0.6218, 0.2269, 0.7705, 0.8326};

	fir_in.add_dist(1);
	filter.init(&fir_in, coeffs);

	// Generate Basis Polynomials
	basis_poly.generate_polys(2);

	//// Set up Simulator
	Simulator sim = Simulator();
	sim.add_basis_poly_set(basis_poly);

	// Add nodes to simulator
	sim.add_node(fir_in);
	sim.add_node(filter);

	sim.set_output_node(filter);
	
	// Set bitwidths
	sim.set_bitwidth(filter, 1);
	sim.set_bitwidth(fir_in, 1);
	
	// Run Sim
	sim.set_sim_params(PCE, tot_sim_steps);
	sim.run_sim(&fir_in);
	// sim.run_sim_anneal(&pll_in, 20.0, 1000);

	filter.tail->print_pwr();
	// out.tail->print_pwr();


	return 0;
}