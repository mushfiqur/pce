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
	int tot_sim_steps = 100;
	
	// Define Netlist
	input_node pll_in(&basis_poly, "pll-in");
	const_node kp(&basis_poly, "kp");
	const_node ki(&basis_poly, "ki");
	const_node k0(&basis_poly, "k0");
	const_node trig_step(&basis_poly, "trig_step");

	//// Phase Detector
	mult_node e_d(&basis_poly, "e_d");
	
	//// PI Filter
	mult_node kp_ed(&basis_poly, "kp_ed");
	mult_node ki_ed(&basis_poly, "ki_ed");
	add_node accum(&basis_poly, "accum");
	delay_node accum_delay(&basis_poly, "accum_delay");
	add_node e_f(&basis_poly, "e_f");

	//// NCO Phase Estimate
	// mult_node k0_ef(&basis_poly, "k0_ef");
	add_node phase_est(&basis_poly, "phase_est");
	add_node phase_est_w_step(&basis_poly, "phase_est_w_step");
	delay_node phase_est_delay(&basis_poly, "phase_est_delay");
	add_node trig_arg(&basis_poly, "trigarg");

	//// NCO
	const_node inverter(&basis_poly, "inverter");
	sine_node sine(&basis_poly, "sine");
	cosine_node cosine(&basis_poly, "cosine");
	mult_node inverted_sine(&basis_poly, "inv_sine");
	delay_node sine_delay(&basis_poly, "sine_delay");

	// Connect/Init Netlist
	pll_in.add_signal(SINE, 1.0/15.0);
	
	pll_in.add_dist(1);
	double snr = 3.0;
	double pwr = 0.5 / std::pow(10.0, snr/10);
	double end = sqrt(3.0*pwr);
	pll_in.set_range(-1.0*end, 1.0*end);

	kp.init(0.2667);
	ki.init(0.0178);
	trig_step.init(2.0*M_PI*(1.0/15.0));

	e_d.init(&pll_in, &inverted_sine);
	kp_ed.init(&e_d, &kp);
	ki_ed.init(&e_d, &ki);
	accum.init(&ki_ed, &accum_delay);
	accum_delay.init(&accum);
	e_f.init(&kp_ed, &accum);

	phase_est.init(&e_f, &phase_est_delay);
	trig_arg.init(&phase_est, &trig_step);
	phase_est_delay.init(&trig_arg);

	inverter.init(-1.0);
	sine.init(&trig_arg);
	cosine.init(&trig_arg);
	sine_delay.init(&sine);
	inverted_sine.init(&sine_delay, &inverter);

	// Generate Basis Polynomials
	basis_poly.generate_polys(2);

	//// Set up Simulator
	Simulator sim = Simulator();
	sim.add_basis_poly_set(basis_poly);

	// Add nodes to simulator
	sim.add_node(pll_in);
	sim.add_node(kp);
	sim.add_node(ki);
	sim.add_node(k0);
	sim.add_node(trig_step);
	sim.add_node(e_d);
	sim.add_node(kp_ed);
	sim.add_node(ki_ed);
	sim.add_node(accum);
	sim.add_node(accum_delay);
	sim.add_node(e_f);
	sim.add_node(phase_est);
	sim.add_node(phase_est_delay);
	sim.add_node(trig_arg);
	sim.add_node(sine);
	sim.add_node(cosine);
	sim.add_node(inverter);
	sim.add_node(inverted_sine);
	sim.add_node(sine_delay);

	sim.set_output_node(cosine);
	
	// Set bitwidths
	sim.set_bitwidth(pll_in, 64);
	sim.set_bitwidth(e_d,    64);
	sim.set_bitwidth(kp_ed,  64);
	sim.set_bitwidth(ki_ed,  64);
	
	// Run Sim
	sim.set_sim_params(PCE, tot_sim_steps);
	// sim.run_sim(&pll_in);
	sim.run_sim_anneal(&pll_in, 20.0, 1000);

	sim.print();

	// sim.add_plot_node(pll_in);
	sim.add_plot_node(cosine);

	sim.plot();

	return 0;
}