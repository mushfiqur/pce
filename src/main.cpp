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
	int correlation = 0;

	////---------------------------------Define Netlist---------------------------------
	input_node input_I(&basis_poly, "I-input");
	input_node input_Q(&basis_poly, "Q-input");

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
	
	////---------------------------------Connect/Init Netlist---------------------------------
	input_I.add_dist(3);
	// input_I.set_range(-1.0*std::sqrt(1.5), 1.0*std::sqrt(1.5));
	input_Q.add_dist(3);
	// input_Q.set_range(-1.0*std::sqrt(1.5), 1.0*std::sqrt(1.5));

	rf_filter_I.init(&input_I, "/home/mushf/pce/filters/rf_filt.txt", correlation);
	rf_filter_Q.init(&input_Q, "/home/mushf/pce/filters/rf_filt.txt", correlation);

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

	// SET OUTPUT NODE
	sim.set_output_node(demod_output);
	
	//// Set bitwidths
	// sim.set_bitwidth(input_I, 24);
	// sim.set_bitwidth(input_Q, 24);
	// // sim.set_bitwidth(a, 8);
	// // sim.set_bitwidth(b, 8);
	// // sim.set_bitwidth(demod_output, 8);
	// // sim.set_bitwidth(rds_channel_squared, 8);
	// // sim.set_bitwidth(e_d, 8);
	// // sim.set_bitwidth(kp_ed, 8);
	// // sim.set_bitwidth(ki_ed, 8);
	// // sim.set_bitwidth(rds_mixed, 8);
	// sim.set_bitwidth(rf_filter_I, 24);
	// sim.set_bitwidth(rf_filter_Q, 24);
	// sim.set_bitwidth(rds_channel_filter, 24);
	// sim.set_bitwidth(rds_carrier_filter, 24);
	// sim.set_bitwidth(rds_baseband_filt, 24);
	// sim.set_bitwidth(rds_rrc_filt, 24);
	
	//// Run Sim
	sim.set_sim_params(PCE, tot_sim_steps);
	sim.run_sim(&input_I);

	// cosine.print_pwr();
	// rds_carrier_filter.print_pwr();
	// sim.add_plot_node(cosine);
	// sim.add_plot_node(rds_rrc_filt);
	// sim.plot();

	rf_filter_I.print(true);
	// demod_output.print(true);
	// demod_output.print_pwr(true);
	// rds_channel_squared.print_pwr(true);
	// rf_filter_I.print_pwr(true);
	// rds_rrc_filt.print_pwr(true);

	return 0;
}