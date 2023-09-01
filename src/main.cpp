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
	int tot_sim_steps = 2;
	
	////---------------------------------Define Netlist---------------------------------
	input_node input_I(&basis_poly, "I_{input}");
	input_node input_Q(&basis_poly, "Q_{input}");

	fir_node rf_filter_I(&basis_poly, "rf-filter-I");
	fir_node rf_filter_Q(&basis_poly, "rf-filter-Q");

	// FM Demodulator //
	const_node c1(&basis_poly, "c1");
	delay_node I_d_1(&basis_poly, "I_d_1");
	delay_node I_d_2(&basis_poly, "I_d_2");
	delay_node Q_d_1(&basis_poly, "Q_d_1");
	delay_node Q_d_2(&basis_poly, "Q_d_2");
	sub_node I_bar(&basis_poly, "i_bar");
	sub_node Q_bar(&basis_poly, "i_bar");
	mult_node a(&basis_poly, "a");
	mult_node b(&basis_poly, "b");
	sub_node pre_scale(&basis_poly, "pre_scale");
	mult_node demod_output(&basis_poly, "demod-out");
	
	// RDS CHANNEL FILTER //
	fir_node rds_channel_filter(&basis_poly, "rds-channel-filt");

	// SQUARING //
	mult_node rds_channel_squared(&basis_poly, "rds-channel-squared");

	// RDS CARRIER FILTER //
	fir_node rds_carrier_filter(&basis_poly, "rds-carrier-filter");

	// PLL //
	input_node pll_signal(&basis_poly, "pll-signal");
	add_node pll_in(&basis_poly, "pll-in");
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

	// MIXER //
	mult_node rds_mixed(&basis_poly, "rds-mixed");

	////---------------------------------Connect/Init Netlist---------------------------------
	input_I.add_dist(3);
	input_Q.add_dist(3);

	rf_filter_I.init(&input_I, "/home/mushf/pce/filters/rf_filt.txt", 0);
	rf_filter_Q.init(&input_Q, "/home/mushf/pce/filters/rf_filt.txt", 0);

	// FM DEMODULATOR //
	c1.init(0.5);
	I_d_1.init(&rf_filter_I);
	I_d_2.init(&I_d_1);
	Q_d_1.init(&rf_filter_Q);
	Q_d_2.init(&Q_d_1);
	I_bar.init(&rf_filter_I, &I_d_2);
	Q_bar.init(&rf_filter_Q, &Q_d_2);
	a.init(&I_d_1, &Q_bar);
	b.init(&Q_d_1, &I_bar);
	pre_scale.init(&a, &b);
	demod_output.init(&pre_scale, &c1);

	// RDS CHANNEL FILTER //
	rds_channel_filter.init(&demod_output, "/home/mushf/pce/filters/rds_channel_filt.txt", 0.0);

	// SQUARING //
	rds_channel_squared.init(&rds_channel_filter, &rds_channel_filter);

	// RDS CARRIER FILTER //
	rds_carrier_filter.init(&rds_channel_squared, "/home/mushf/pce/filters/rds_carrier_filt.txt", 0.0);

	// PLL //
	pll_signal.add_signal(SINE, 1.0/15.0);
	pll_in.init(&rds_carrier_filter, &pll_signal);
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

	// RDS MIXED //
	// const_node c_debug(&basis_poly, "c-debug");
	// c_debug.init(1.0);
	// rds_mixed.init(&c_debug, &cosine);
	rds_mixed.init(&rds_channel_filter, &cosine);

	//// Generate Basis Polynomials
	basis_poly.generate_polys(2);

	////---------------------------------Set up Simulator-----------------------
	Simulator sim = Simulator();
	sim.add_basis_poly_set(basis_poly);

	// Add nodes to simulator
	sim.add_node(input_I);
	sim.add_node(input_Q);
	sim.add_node(rf_filter_I);
	sim.add_node(rf_filter_Q);
	sim.add_node(c1);
	sim.add_node(I_d_1);
	sim.add_node(I_d_2);
	sim.add_node(Q_d_1);
	sim.add_node(Q_d_2);
	sim.add_node(I_bar);
	sim.add_node(Q_bar);
	sim.add_node(a);
	sim.add_node(b);
	sim.add_node(pre_scale);
	sim.add_node(demod_output);
	sim.add_node(rds_channel_filter);
	sim.add_node(rds_channel_squared);
	sim.add_node(rds_carrier_filter);
	sim.add_node(pll_signal);
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
	sim.add_node(rds_mixed);
	// sim.add_node(c_debug);
	
	// SET OUTPUT NODE
	sim.set_output_node(cosine);
	
	//// Set bitwidths
	// sim.set_bitwidth(rf_filter, 1);
	
	//// Run Sim
	sim.set_sim_params(PCE, tot_sim_steps);
	sim.run_sim(&input_I);

	// cosine.print_pwr();
	// rds_carrier_filter.print_pwr();
	sim.add_plot_node(cosine);
	sim.add_plot_node(rds_mixed);
	sim.plot();



	return 0;
}