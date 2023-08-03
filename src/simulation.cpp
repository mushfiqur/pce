#include "../include/simulation.h"

Simulator::Simulator(){
	this->pce_sim_done = false;
	this->mc_sim_done = false;

	std::random_device rd;
	this->mt = std::mt19937(rd());
	this->real_dist = std::uniform_real_distribution<double>(0.0, 1.0);
}

void Simulator::set_sim_params(SimType sim_t, int tot_sim_steps, int mc_samples){
	this->tot_sim_steps = tot_sim_steps;
	this->mc_samples = mc_samples;
	this->sim_t = sim_t;
}

void Simulator::add_node(dfg_node& n){
	this->nodes_arr.push_back(&n);
	// if(n.t == INPUT || n.t == MULT){
	// 	this->curr_config.push_back({.n = &n, .bitwidth=-1});
	// }
}

void Simulator::set_node_sim_params(){
	for(int i = 0; i < this->nodes_arr.size(); i++){
		this->nodes_arr[i]->set_sim_params(this->tot_sim_steps, this->mc_samples, this->bp_set_ptr->set_size, this->sim_t);
	}
}

void Simulator::set_bitwidth(dfg_node& n, int bitwidth){
	for(int i = 0; i < this->curr_solution.size(); i++){
		if(this->curr_solution[i].n == &n){
			this->curr_solution[i].bitwidth == bitwidth;
			return;
		}
	}
	this->curr_solution.push_back({.n = &n, .bitwidth=bitwidth});
	// std::clog << "ERROR: cannot set bitwidth for node " << n.label << std::endl;
}

void Simulator::set_output_node(dfg_node& n){
	this->output_node = &n;
}

void Simulator::calc_bitwidths(){
	// TODO: If node is part of a circular graph, then bitwidth is INF
	//        else, node bitwidth is function of lhs and rhs
}

void Simulator::initialize(dfg_node* n){
	this->head = n;
	this->int_dist = std::uniform_int_distribution<int>(0, curr_solution.size() - 1);
	
	this->set_node_sim_params();

	if(this->sim_t == MONTE_CARLO){
		this->mc_sim_done = true;
	}
	else{
		this->pce_sim_done = true;
	}

	// Propagate signal
	this->propagate_coeffs();

	// Store signal coeffs
	for(int i = 0; i < this->nodes_arr.size(); i++){
		this->nodes_arr[i]->save_signal_polys();
	}

	// Set bitwidths
	for(int i = 0; i < this->curr_solution.size(); i++){
		this->curr_solution[i].n->set_bitwidth(this->curr_solution[i].bitwidth);
	}

	//// Regenerate basis_polys
	bp_set_ptr->regenerate_polys(3);
	
	//// Reorder signal coeffs in new basis_poly indexing
	for(int i = 0; i < nodes_arr.size(); i++){
		nodes_arr[i]->reorder_signal_polys();
	}

	// Propagate signal + noise
	this->set_node_sim_params();
	this->propagate_coeffs();

	// Remove signal component
	for(int n = 0; n < this->nodes_arr.size(); n++){
		this->nodes_arr[n]->remove_signal_component();
	}

	// Save power
	this->curr_sol_noise_pwr = this->output_node->get_pwr();
}

double Simulator::try_solution(std::vector<bitwidth_config>& proposed_sol){
	// Set bitwidths
	for(int i = 0; i < proposed_sol.size(); i++){
		proposed_sol[i].n->set_bitwidth(proposed_sol[i].bitwidth);
	}

	// Propagate signal + noise
	this->set_node_sim_params();
	this->propagate_coeffs();

	// Remove signal component
	for(int n = 0; n < this->nodes_arr.size(); n++){
		this->nodes_arr[n]->remove_signal_component();
	}

	// Save noise power
	return this->output_node->get_pwr();
}

void Simulator::print(){
	std::cout << "\n(" << this->curr_sol_noise_pwr << ") { \n";
	for(int i = 0; i < this->curr_solution.size() - 1; i++){
		std::cout << "\t" << this->curr_solution[i].n->label << ": " << this->curr_solution[i].bitwidth << ", \n";
	}
	std::cout << "\t" << this->curr_solution[this->curr_solution.size() - 1].n->label << ": " << this->curr_solution[this->curr_solution.size() - 1].bitwidth << std::endl;
	std::cout << "}\n";
}

void Simulator::run_sim_anneal(dfg_node* n, double sig_pwr, double tgt_snr, int tot_iters){
		// Vars
	std::vector<bitwidth_config> proposed_sol;
	double prop_sol_noise_pwr;

	double tgt_noise_pwr = sig_pwr / std::pow(10.0, tgt_snr / 10.0);
	double temperature = 0.1;
	
	// Load initial solution
	this->initialize(n);
	
	std::cout << "Tgt: " << tgt_noise_pwr << std::endl;
	std::cout << this->curr_sol_noise_pwr << std::endl;

	for(int iter = 0; iter < tot_iters; iter++){
		// Get neighboring solution
		proposed_sol = get_neighbour();

		// Try solution
		prop_sol_noise_pwr = try_solution(proposed_sol);

		// Evaluate solution
		if(prop_sol_noise_pwr <= tgt_noise_pwr){
			if(prop_sol_noise_pwr > curr_sol_noise_pwr){
				this->curr_solution = proposed_sol;
				curr_sol_noise_pwr = prop_sol_noise_pwr;
				// std::cout << "[" << iter << "] " << "{" << temperature << "} ";
				// this->print();
			}
			else{
				// std::cout << "Probability of accept: " << std::exp((prop_sol_noise_pwr - curr_sol_noise_pwr) / temperature) << std::endl;
				if (std::exp((prop_sol_noise_pwr - curr_sol_noise_pwr) / temperature) >= this->real_dist(this->mt) ){
					this->curr_solution = proposed_sol;
					curr_sol_noise_pwr = prop_sol_noise_pwr;
					// std::cout << "LUCKY ";
					// std::cout << "[" << iter << "] " << "{" << temperature << "} ";
					// this->print();
				}
			}
		}

		// Update temperature
		// temperature = 1.0 / std::log10(2 + iter);
		temperature = 0.01 * temperature;
	}

	// print();
}

void Simulator::run_sim(dfg_node* n){
	this->head = n;
	this->int_dist = std::uniform_int_distribution<int>(0, curr_solution.size() - 1);
	
	this->set_node_sim_params();

	if(this->sim_t == MONTE_CARLO){
		this->mc_sim_done = true;
	}
	else{
		this->pce_sim_done = true;
	}

	// Propagate signal
	this->propagate_coeffs();

	if(this->curr_solution.size() != 0){
		// Store signal coeffs
		for(int i = 0; i < this->nodes_arr.size(); i++){
			this->nodes_arr[i]->save_signal_polys();
		}

		// Set bitwidths
		for(int i = 0; i < this->curr_solution.size(); i++){
			this->curr_solution[i].n->set_bitwidth(this->curr_solution[i].bitwidth);
		}

		//// Regenerate basis_polys
		bp_set_ptr->regenerate_polys(2);
		
		//// Reorder signal coeffs in new basis_poly indexing
		for(int i = 0; i < nodes_arr.size(); i++){
			nodes_arr[i]->reorder_signal_polys();
		}

		// Propagate signal + noise
		this->set_node_sim_params();
		this->propagate_coeffs();

		// Remove signal component
		for(int n = 0; n < this->nodes_arr.size(); n++){
			this->nodes_arr[n]->remove_signal_component();
		}

		// Save power
		this->curr_sol_noise_pwr = this->output_node->get_pwr();
	}
}

std::vector<bitwidth_config> Simulator::get_neighbour(){
	std::vector<bitwidth_config> proposed_sol(this->curr_solution);
	int rand_int = this->int_dist(this->mt);

	
	while(proposed_sol[rand_int].bitwidth == 1){
		rand_int = this->int_dist(this->mt);
	}
	proposed_sol[rand_int].bitwidth -= 1;

    // std::srand((unsigned)time(0)); 
	// if(proposed_sol[rand_int].bitwidth > 1){
	// 	proposed_sol[rand_int].bitwidth -= 1;;
	// }
	// proposed_sol[rand_int].bitwidth = (std::rand()%32)+1;;

	return proposed_sol;
}

void Simulator::propagate_coeffs(){
	std::deque<dfg_node*> q;

	for(int curr_timestamp = 0; curr_timestamp < this->tot_sim_steps; curr_timestamp++){
		q.push_back(this->head);

		dfg_node* curr_node;
		while(q.size() != 0){
			curr_node = q.front();

			if(curr_node->t == CONST || curr_node->t == DELAY || curr_node->t == INPUT_SIGNAL || curr_node->t == INPUT_NOISE || curr_node->node_args_ready(curr_timestamp)){
				// std::cout << "[" << curr_timestamp << "] Processing " << curr_node->label << std::endl;
				curr_node->process(curr_timestamp);
				q.pop_front();
				for(int i = 0; i < curr_node->next_nodes.size(); i++){
					if( curr_node->next_nodes[i]->last_exec_time < curr_timestamp && std::find(q.begin(), q.end(), curr_node->next_nodes[i]) == q.end()){
						q.push_back(curr_node->next_nodes[i]);
					}
				}
			}
			else{
				if(curr_node->lhs != nullptr && curr_node->lhs->last_exec_time < curr_timestamp){
					if(std::find(q.begin(), q.end(), curr_node->lhs) == q.end()){
						q.push_front(curr_node->lhs);
					}
					else{
						q.erase(std::find(q.begin(), q.end(), curr_node->lhs));
						q.push_front(curr_node->lhs);
					}
				}
				
				if(curr_node->rhs != nullptr && curr_node->rhs->last_exec_time < curr_timestamp){
					if(std::find(q.begin(), q.end(), curr_node->rhs) == q.end()){
						q.push_front(curr_node->rhs);
					}
					else{
						q.erase(std::find(q.begin(), q.end(), curr_node->rhs));
						q.push_front(curr_node->rhs);
					}
				}

				// if(curr_node->rhs != nullptr &&
				// curr_node->rhs->last_exec_time < curr_timestamp &&
				// std::find(q.begin(), q.end(), curr_node->rhs) == q.end())
				// {
				// 	q.push_front(curr_node->rhs);
				// }
			}
		}

		// std::cout << "---------" << std::endl;
	}
}

void Simulator::add_plot_node(dfg_node& n){
	this->nodes_to_plot.push_back(&n);
}

void Simulator::plot(){
	std::vector<double> mean_arr(this->tot_sim_steps, 0.0);
	std::vector<double> var_arr(this->tot_sim_steps, 0.0);

	std::vector<std::vector<double>> plot_data_arr;
	std::string plot_cmd;

	double mc_var = 0.0;
	double pce_var = 0.0;
	double sim_err = 0.0;

	int num_max_per_row = 2;
	int nrows, ncols;
	if(this->nodes_to_plot.size() < num_max_per_row){
		nrows = this->nodes_to_plot.size();
		ncols = 1;
	}
	else{
		nrows = num_max_per_row;
		ncols = (int)std::ceil(this->nodes_to_plot.size() / (float)num_max_per_row);
	}

	if(nrows == 0){
		nrows = 1;
	}
	if(ncols == 0){
		ncols = 1;
	}

	// Plotting
	Gnuplot gp("gnuplot");    
	gp << "set terminal qt size 1000,600 title 'Simulation Results'\n";
	// gp << "set terminal svg size 1000,600\n";
	gp << "set output '/home/mushf/pce/plots/results.svg'\n";
	gp << "set multiplot layout " << std::to_string(nrows) << "," << std::to_string(ncols) << " rowsfirst\n";

	for(int i = 0; i < this->nodes_to_plot.size(); i++){
		plot_cmd = "plot ";
		plot_data_arr.clear();
		
		gp << "set title 'Node " << this->nodes_to_plot[i]->label.c_str() << " Simulation Results' font ',18' \n";
		
		if(this->mc_sim_done){
			this->nodes_to_plot[i]->get_mc_stats(mean_arr, var_arr);
			mc_var = var_arr[this->tot_sim_steps - 1];
			// plot_cmd.append(" '-' with linespoints dt \".\" lc rgb 'blue'");
			// plot_cmd.append(" title 'MC Mean',");
			// plot_data_arr.push_back(mean_arr);
			
			plot_cmd.append(" '-' with linespoints dt \".\" lc rgb 'blue'");
			plot_cmd.append(" title 'MC Var.',");
			plot_data_arr.push_back(var_arr);
		}
		if(this->pce_sim_done){
			this->nodes_to_plot[i]->get_pce_stats(mean_arr, var_arr, *bp_set_ptr);
			pce_var = var_arr[this->tot_sim_steps - 1];

			// plot_cmd.append(" '-' with lines lc rgb 'red'");
			// plot_cmd.append(" title 'PCE Mean',");
			// plot_data_arr.push_back(mean_arr);
			
			plot_cmd.append(" '-' with lines lc rgb '#005A32'");
			plot_cmd.append(" title 'PCE Var.',");
			plot_data_arr.push_back(var_arr);
		}
		gp << "set key outside\n";
		gp << "set grid\n";

		if(plot_cmd.back() == ','){
			plot_cmd.pop_back();
		}
		plot_cmd.append("\n");

		gp << "set xlabel 'Iteration' font ',16'\n";
		gp << "set ylabel 'Moment' font ',16'\n";

		gp << plot_cmd;

		for(int i = 0; i < plot_data_arr.size(); i++){
			gp.send(plot_data_arr[i]);
		}

		if(this->pce_sim_done && this->mc_sim_done){
			sim_err = (std::abs(pce_var - mc_var) / mc_var) * 100.0;
			std::clog << "[" << this->nodes_to_plot[i]->label << "]" << std::endl;
			std::clog << "PCE Variance: " << pce_var << std::endl;
			std::clog << "MC Variance: " << mc_var << std::endl;
			std::clog << "Simulation Error: " << sim_err << "%\n" << std::endl;
		}
	}
	

	gp << "unset multiplot\n";

	gp << "pause mouse close\n";

}

void Simulator::add_basis_poly_set(BasisPolySet& bp_set){
	this->bp_set_ptr = &bp_set;
}