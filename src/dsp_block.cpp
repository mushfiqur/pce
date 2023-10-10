#include "../include/dsp_block.h"

sine_node::sine_node(BasisPolySet* bp_set, std::string label) : dfg_node(SINE_BLOCK, bp_set, label) { 

}

void sine_node::init(dfg_node* arg_node){
	this->head->lhs = arg_node->tail;
	this->head->rhs = nullptr;

	arg_node->tail->add_next_node(this);
}

void sine_node::process(int curr_timestamp){
	this->last_exec_time = curr_timestamp;

	// d = arg_node - mean(arg_node)
	std::vector<double> d = this->head->lhs->pce_coeffs[curr_timestamp];
	double center = d[0];
	d[0] = 0.0;

	// Top_a = d*d
	std::vector<double> top_a = std::vector<double>(d.size(), 0.0);
	for(int k = 0; k < top_a.size(); k++){
		for(int i = 0; i < d.size(); i++){
			for(int j = 0; j < d.size(); j++){
				top_a[k] += d[i] * d[j] * this->bp_set_ptr->expt_table[i][j][k];
			}
		}
		top_a[k] /= this->bp_set_ptr->poly_sqr_expt[k];
	}

	// Top_b = top_a*d
	std::vector<double> top_b = std::vector<double>(d.size(), 0.0);
	for(int k = 0; k < top_b.size(); k++){
		for(int i = 0; i < top_a.size(); i++){
			for(int j = 0; j < d.size(); j++){
				top_b[k] += top_a[i] * d[j] * this->bp_set_ptr->expt_table[i][j][k];
			}
		}
		top_b[k] /= this->bp_set_ptr->poly_sqr_expt[k];
	}
	
	// Top_b = top_b * top_const
	double top_const = (-1.0/6.0)*std::cos(center);
	for(int i = 0; i < top_b.size(); i++){
		top_b[i] = top_b[i] * top_const;
	}

	// Mid_a = top_a * mid_const
	std::vector<double> mid_a = std::vector<double>(d.size(), 0.0);
	double mid_const = (-1.0/2.0)*std::sin(center);
	for(int i = 0; i < mid_a.size(); i++){
		mid_a[i] = top_a[i] * mid_const;
	}

	// Bot_a = d * bot_const
	double bot_const = std::cos(center);
	std::vector<double> bot_a = std::vector<double>(d.size(), 0.0);
	for(int i = 0; i < bot_a.size(); i++){
		bot_a[i] = d[i] * bot_const;
	}

	// out = end_const + sum(top_b, mid_a, bot_a)
	double end_const = std::sin(center);
	for(int i = 0; i < this->pce_coeffs[curr_timestamp].size(); i++){
		this->pce_coeffs[curr_timestamp][i] = (top_b[i] + mid_a[i] + bot_a[i]);
	}
	this->pce_coeffs[curr_timestamp][0] += end_const;
}

void sine_node::set_sim_params(int tot_sim_steps, int mc_samples, int basis_set_size, SimType sim_type){
	this->last_exec_time = -1;
	this->sim_type = sim_type;

	this->pce_coeffs = std::vector<std::vector<double>>(tot_sim_steps, std::vector<double>(basis_set_size, 0.0));
}

void sine_node::print(bool print_last){
	this->tail->print(print_last);
}

void sine_node::set_bitwidth(int width){
	// Not supported
	return;
}

void sine_node::reorder_signal_polys(){
	dfg_node::reorder_signal_polys();
}

void sine_node::save_signal_polys(){
	dfg_node::save_signal_polys();
}

void sine_node::remove_signal_component(){
	dfg_node::remove_signal_component();
}

double sine_node::get_pwr(){
	double alpha = 0.05;
	double alpha2 = 0.1;

	dfg_node* n = this->tail;

	std::vector<double> out(n->pce_coeffs.size(), 0.0);
	std::vector<double> arr(n->pce_coeffs.size(), 0.0);

	for(int i = 0; i < arr.size(); i++){
		arr[i] = n->pce_coeffs[t][i] * n->pce_coeffs[t][i] * this->bp_set_ptr->poly_sqr_expt[i];
	}

	for(int n = 0; n < out.size(); n++){
		if(n - 1 > 0){
			out[n] = (1.0 - alpha)*out[n-1] + alpha*arr[n];
		}
		else{
			out[n] = alpha*arr[n];
		}
	}

	for(int n = 0; n < out.size(); n++){
		if(n - 1 > 0){
			out[n] = (1.0 - alpha2)*out[n-1] + alpha2*out[n];
		}
		else{
			out[n] = alpha2*out[n];
		}
	}

	double max = 0.0;
	for(int i = 0; i < out.size(); i++){
		if(out[i] > max){
			max = out[i];
		}
	}

	return max;
	// return out[out.size() - 1];
}

////////////////////////////////////////////////////////////////////////////////////

cosine_node::cosine_node(BasisPolySet* bp_set, std::string label) : dfg_node(COSINE_BLOCK, bp_set, label) { 

}

void cosine_node::init(dfg_node* arg_node){
	this->head->lhs = arg_node->tail;
	this->head->rhs = nullptr;

	arg_node->tail->add_next_node(this);
}

void cosine_node::process(int curr_timestamp){
	
	this->last_exec_time = curr_timestamp;

	// d = arg_node - mean(arg_node)
	std::vector<double> d = this->head->lhs->pce_coeffs[curr_timestamp];
	double center = d[0];
	d[0] = 0.0;

	// Top_a = d*d
	std::vector<double> top_a = std::vector<double>(d.size(), 0.0);
	for(int k = 0; k < top_a.size(); k++){
		for(int i = 0; i < d.size(); i++){
			for(int j = 0; j < d.size(); j++){
				top_a[k] += d[i] * d[j] * this->bp_set_ptr->expt_table[i][j][k];
			}
		}
		top_a[k] /= this->bp_set_ptr->poly_sqr_expt[k];
	}

	// Top_b = top_a*d
	std::vector<double> top_b = std::vector<double>(d.size(), 0.0);
	for(int k = 0; k < top_b.size(); k++){
		for(int i = 0; i < top_a.size(); i++){
			for(int j = 0; j < d.size(); j++){
				top_b[k] += top_a[i] * d[j] * this->bp_set_ptr->expt_table[i][j][k];
			}
		}
		top_b[k] /= this->bp_set_ptr->poly_sqr_expt[k];
	}
	
	// Top_b = top_b * top_const
	double top_const = (1.0/6.0)*std::sin(center);
	for(int i = 0; i < top_b.size(); i++){
		top_b[i] = top_b[i] * top_const;
	}

	// Mid_a = top_a * mid_const
	std::vector<double> mid_a = std::vector<double>(d.size(), 0.0);
	double mid_const = (-1.0/2.0)*std::cos(center);
	for(int i = 0; i < mid_a.size(); i++){
		mid_a[i] = top_a[i] * mid_const;
	}

	// Bot_a = d * bot_const
	double bot_const = -1.0*std::sin(center);
	std::vector<double> bot_a = std::vector<double>(d.size(), 0.0);
	for(int i = 0; i < bot_a.size(); i++){
		bot_a[i] = d[i] * bot_const;
	}

	// out = end_const + sum(top_b, mid_a, bot_a)
	double end_const = std::cos(center);
	for(int i = 0; i < this->pce_coeffs[curr_timestamp].size(); i++){
		this->pce_coeffs[curr_timestamp][i] = (top_b[i] + mid_a[i] + bot_a[i]);
	}
	this->pce_coeffs[curr_timestamp][0] += end_const;
}

void cosine_node::set_sim_params(int tot_sim_steps, int mc_samples, int basis_set_size, SimType sim_type){
	this->last_exec_time = -1;
	this->sim_type = sim_type;

	this->pce_coeffs = std::vector<std::vector<double>>(tot_sim_steps, std::vector<double>(basis_set_size, 0.0));
}

void cosine_node::print(bool print_last){
	this->tail->print(print_last);
}

void cosine_node::set_bitwidth(int width){
	// Not supported
	return;
}

void cosine_node::reorder_signal_polys(){
	dfg_node::reorder_signal_polys();
}

void cosine_node::save_signal_polys(){
	dfg_node::save_signal_polys();
}

void cosine_node::remove_signal_component(){
	dfg_node::remove_signal_component();
}

double cosine_node::get_pwr(){
	double alpha = 0.05;
	double alpha2 = 0.1;

	dfg_node* n = this->tail;

	std::vector<double> out(n->pce_coeffs.size(), 0.0);
	std::vector<double> arr(n->pce_coeffs.size(), 0.0);

	for(int t = 0; t < arr.size(); t++){
		arr[t] = 0.0;
		for(int i = 0; i < n->pce_coeffs[t].size(); i++){
			arr[t] += n->pce_coeffs[t][i] * n->pce_coeffs[t][i] * this->bp_set_ptr->poly_sqr_expt[i];
		}
	}

	for(int n = 0; n < out.size(); n++){
		if(n - 1 > 0){
			out[n] = (1.0 - alpha)*out[n-1] + alpha*arr[n];
		}
		else{
			out[n] = alpha*arr[n];
		}
	}

	for(int n = 0; n < out.size(); n++){
		if(n - 1 > 0){
			out[n] = (1.0 - alpha2)*out[n-1] + alpha2*out[n];
		}
		else{
			out[n] = alpha2*out[n];
		}
	}

	double max = 0.0;
	for(int i = 0; i < out.size(); i++){
		if(out[i] > max){
			max = out[i];
		}
	}

	return max;
	// return out[out.size() - 1];
}

void cosine_node::get_pce_stats(std::vector<double>& mean_arr, std::vector<double>& var_arr){
	double alpha = 0.05;
	double alpha2 = 0.1;

	dfg_node* n = this->tail;

	for(int t = 0; t < var_arr.size(); t++){
		var_arr[t] = 0.0;
		mean_arr[t] = 0.0;
		for(int i = 0; i < n->pce_coeffs[t].size(); i++){
			var_arr[t] += n->pce_coeffs[t][i] * n->pce_coeffs[t][i] * this->bp_set_ptr->poly_sqr_expt[i];
			mean_arr[t] += n->pce_coeffs[t][i] * this->bp_set_ptr->poly_expt[i];
		}
	}

	for(int n = 0; n < var_arr.size(); n++){
		if(n - 1 > 0){
			var_arr[n] = (1.0 - alpha)*var_arr[n-1] + alpha*var_arr[n];
		}
		else{
			var_arr[n] = alpha*var_arr[n];
		}
	}

	for(int n = 0; n < var_arr.size(); n++){
		if(n - 1 > 0){
			var_arr[n] = (1.0 - alpha2)*var_arr[n-1] + alpha2*var_arr[n];
		}
		else{
			var_arr[n] = alpha2*var_arr[n];
		}
	}

}

////////////////////////////////////////////////////////////////////////////////////

collate_node::collate_node(BasisPolySet* bp_set, std::string label) : dfg_node(DELAY, bp_set, label){
	var* v = new var(this->bp_set_ptr->get_new_var_id());
	v->init_uniform_var(-1, 1);
	this->v = v;
	this->bp_set_ptr->add_variable(v);
}

void collate_node::init(dfg_node* arg_node){
	this->head->lhs = arg_node->tail;
	this->head->rhs = nullptr;

	arg_node->tail->add_next_node(this);
}

void collate_node::set_bitwidth(int width){
	// Not supported
	return;
}

void collate_node::set_sim_params(int tot_sim_steps, int mc_samples, int basis_set_size, SimType sim_type){
	this->sim_type = sim_type;
	this->last_exec_time = -1;

	this->pce_coeffs = std::vector<std::vector<double>>(tot_sim_steps, std::vector<double>(basis_set_size, 0.0));
}

void collate_node::process(int curr_timestamp){
	this->last_exec_time = curr_timestamp;

	// Get mean and var from prev timestamp
	double mean = 0.0;
	double var = 0.0;

	if(curr_timestamp - 1 >= 0){
		mean = lhs->pce_coeffs[curr_timestamp - 1][0];

		for(int i = 0; i < this->lhs->pce_coeffs[curr_timestamp - 1].size(); i++){
			var  += lhs->pce_coeffs[curr_timestamp - 1][i] * lhs->pce_coeffs[curr_timestamp - 1][i] * this->bp_set_ptr->poly_sqr_expt[i];
		}

		this->pce_coeffs[curr_timestamp][0] = mean;
		this->pce_coeffs[curr_timestamp][this->bp_set_ptr->get_var_idx(this->v->id)] = std::sqrt(3.0 * (var - (mean * mean)));
	}
	else{
		for(int i = 0; i < this->pce_coeffs[curr_timestamp].size(); i++){
			this->pce_coeffs[curr_timestamp][i] = 0.0;
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////

fir_node::fir_node(BasisPolySet* bp_set, std::string label) : dfg_node(FIR_BLOCK, bp_set, label){
	// var* v = new var(1, -1, this->bp_set_ptr->get_new_var_id());
	// this->v = v;
	// this->bp_set_ptr->add_variable(v);

	this->sum_coeffs = 0.0;
	this->sum_sqr_coeffs = 0.0;
}

void fir_node::init(dfg_node* arg_node, std::string filename, int correlation_dist){
	this->head->lhs = arg_node->tail;
	this->head->rhs = nullptr;

	arg_node->tail->add_next_node(this);

	std::ifstream f_coeffs(filename, std::ios::in);
	double num_read = 0.0;
	
	while(f_coeffs >> num_read){
		this->filt_coeffs.push_back(num_read);
		this->sum_coeffs += num_read;
		this->sum_sqr_coeffs += num_read * num_read;
	}
	
	// for(int i = 0; i < filt_coeffs.size(); i++){
	// 	for(int j = 0; j < filt_coeffs.size(); j++){
	// 		if(std::abs(i - j) <= correlation_dist){
	// 			this->sum_sqr_coeffs += std::abs(filt_coeffs[i] * filt_coeffs[j]);
	// 		}
	// 		// else{
	// 		// 	this->sum_sqr_coeffs += correlation * filt_coeffs[i] * filt_coeffs[j];
	// 		// }
			
	// 	}
	// }
}

void fir_node::init(dfg_node* arg_node, std::vector<double>& coeffs){
	this->head->lhs = arg_node->tail;
	this->head->rhs = nullptr;
	this->n_node = nullptr;

	arg_node->tail->add_next_node(this);

	for(int i = 0; i < coeffs.size(); i++){
		this->filt_coeffs.push_back(coeffs[i]);
	}

	for(int i = 0; i < coeffs.size(); i++){

		for(int j = 0; j < coeffs.size(); j++){
			this->sum_sqr_coeffs += coeffs[i] * coeffs[j];
		}
	}
}

void fir_node::set_sim_params(int tot_sim_steps, int mc_samples, int basis_set_size, SimType sim_type){
	this->sim_type = sim_type;
	this->last_exec_time = -1;

	// Initialize pce array for this node
	this->pce_coeffs = std::vector<std::vector<double>>(tot_sim_steps, std::vector<double>(basis_set_size, 0.0));

	
	if(this->tail != this){
		this->tail->set_sim_params(tot_sim_steps, mc_samples, basis_set_size, sim_type);
		this->tail->rhs->set_sim_params(tot_sim_steps, mc_samples, basis_set_size, sim_type);

		// dfg_node* n = this->tail->rhs;	// noise node
		// noise_node* n = static_cast<noise_node*>(this->tail->rhs);
		

		// Initialize pce array for noise summation node
		// this->tail->pce_coeffs = std::vector<std::vector<double>>(tot_sim_steps, std::vector<double>(basis_set_size, 0.0));
		// Initialize pce array for noise node
		// this->n_node->pce_coeffs = std::vector<std::vector<double>>(tot_sim_steps, std::vector<double>(basis_set_size, 0.0));

		
		// for(int t = 0; t < this->n_node->pce_coeffs.size(); t++){
		// 	this->n_node->pce_coeffs[t][this->bp_set_ptr->get_var_idx(this->n_node->vars[t % this->vars.size()]->id)] = (std::sqrt(3.0 * this->filt_coeffs.size()) * std::pow(2.0, -1*bitwidth)) / 6.0;
		// }

	}
}

void fir_node::set_bitwidth(int width){
	this->bitwidth = width;

	if(this->tail == this){
		int n_noise_id = this->bp_set_ptr->get_new_var_id();
		std::string lbl = "n" + std::to_string(n_noise_id);
		noise_node* n_noise = new noise_node(GAUSSIAN, this->vars.size(), this->bp_set_ptr, n_noise_id, lbl);;
		n_noise->bitwidth = width;
		this->n_node = n_noise;

		this->tail = new add_node(this->bp_set_ptr, this->label + "_n");
		this->tail->lhs = this;
		this->tail->rhs = n_noise;

		for(int i = 0; i < this->next_nodes.size(); i++){
			this->tail->add_next_node(this->next_nodes[i]);

			if(this->next_nodes[i]->head->lhs == this){
				this->next_nodes[i]->head->lhs = this->tail;
			}

			if(this->next_nodes[i]->head->rhs == this){
				this->next_nodes[i]->head->rhs = this->tail;
			}

			for(int j = 0; j < this->next_nodes[i]->prev_nodes.size(); j++){
				if(this->next_nodes[i]->prev_nodes[j]->label == this->label){
					this->next_nodes[i]->prev_nodes.erase(this->next_nodes[i]->prev_nodes.begin() + j);
				}
			}
		}
		this->next_nodes.clear();

		this->add_next_node(this->tail);
		n_noise->add_next_node(this->tail);
	}
	else{
		this->tail->rhs->set_bitwidth(width);	
	}
}

void fir_node::print(bool print_last){
	this->tail->print(print_last);
}

void fir_node::add_dist(int num_rand_vars){
	var* v;
	for(int i = 0; i < num_rand_vars; i++){
		v = new var(this->bp_set_ptr->get_new_var_id());
		v->init_gaussian_var(0.0, 1.0);
		this->bp_set_ptr->add_variable(v);
		this->vars.push_back(v);
	}
	
}

// void fir_node::remove_signal_component(){
// 	for(int t = 0; t < this->pce_coeffs.size(); t++){
// 		for(int i = 0; i < this->pce_coeffs[t].size(); i++){
// 			// this->pce_coeffs[t][i] -= this->signal_coeffs[t][i];
// 			this->pce_coeffs[t][i] = std::sqrt(std::pow(this->pce_coeffs[t][i], 2.0) - std::pow(this->signal_coeffs[t][i], 2.0));
// 			if(this->tail != this){
// 				// this->tail->pce_coeffs[t][i] = std::sqrt(std::pow(this->tail->pce_coeffs[t][i], 2.0) - std::pow(this->signal_coeffs[t][i], 2.0));
// 				this->tail->pce_coeffs[t][i] -= this->signal_coeffs[t][i];
// 			}
// 		}
// 	}
// }

void fir_node::process(int curr_timestamp){
	this->last_exec_time = curr_timestamp;
	int noise_rv_ctr;

	double signal_variance = 0.0;		// Add to FIR sig id's coefficient
	double noise_variance = 0.0;		// Add to fir noise id's coefficient
	double mix_variance = 0.0;			// Add to fir product poly's coefficient

	for(int i = 1; i < this->lhs->pce_coeffs[curr_timestamp].size(); i++){
		noise_rv_ctr = 0;

		for(int j = 0; j < this->bp_set_ptr->basis_polys[i]->var_ids_contained.size(); j++){
			for(int k = 0; k < this->bp_set_ptr->noise_vars.size(); k++){
				if(this->bp_set_ptr->basis_polys[i]->var_ids_contained[j] == this->bp_set_ptr->noise_vars[k]->id){
					noise_rv_ctr++;
				}
			}
		}

		if(noise_rv_ctr == 0){
			// purely signal node
			// std::cout << "Poly " << i << " is signal poly" << std::endl;
			signal_variance += std::pow(this->lhs->pce_coeffs[curr_timestamp][i], 2.0) * this->bp_set_ptr->poly_sqr_expt[i];
		}
		else if(noise_rv_ctr == this->bp_set_ptr->basis_polys[i]->var_ids_contained.size()){
			// purely noise node
			// std::cout << "Poly " << i << " is noise poly" << std::endl;
			noise_variance += std::pow(this->lhs->pce_coeffs[curr_timestamp][i], 2.0) * this->bp_set_ptr->poly_sqr_expt[i];
		}
		else{
			// mixed node
			// std::cout << "Poly " << i << " is mixed poly" << std::endl;
			mix_variance += std::pow(this->lhs->pce_coeffs[curr_timestamp][i], 2.0) * this->bp_set_ptr->poly_sqr_expt[i];
		}
	}

	this->pce_coeffs[curr_timestamp][0] = this->sum_coeffs * lhs->pce_coeffs[curr_timestamp][0];
	
	int sig_idx = this->vars[curr_timestamp % this->vars.size()]->id;
	this->pce_coeffs[curr_timestamp][this->bp_set_ptr->get_var_idx(sig_idx)] = std::sqrt(this->sum_sqr_coeffs * signal_variance);

	if(this->n_node != nullptr){
		int noise_idx = this->n_node->vars[curr_timestamp % this->vars.size()]->id;	
		this->n_node->pce_coeffs[curr_timestamp][this->bp_set_ptr->get_var_idx(noise_idx)] = ((std::sqrt(3.0 * this->filt_coeffs.size()) * std::pow(2.0, -1*n_node->bitwidth)) / 6.0) + std::sqrt(this->sum_sqr_coeffs * noise_variance);	
		this->pce_coeffs[curr_timestamp][this->bp_set_ptr->get_var_idx(sig_idx, noise_idx)] = std::sqrt(this->sum_sqr_coeffs * mix_variance);
	}
	return;

	// double var = 0.0;
	// for(int i = 1; i < this->lhs->pce_coeffs[curr_timestamp].size(); i++){
	// 	var += std::pow(this->lhs->pce_coeffs[curr_timestamp][i], 2.0) * this->bp_set_ptr->poly_sqr_expt[i];
	// }

	// var *= this->sum_sqr_coeffs;

	// this->pce_coeffs[curr_timestamp][0] = this->sum_coeffs * lhs->pce_coeffs[curr_timestamp][0];
	// this->pce_coeffs[curr_timestamp][this->bp_set_ptr->get_var_idx(this->vars[curr_timestamp % this->vars.size()]->id)] = std::sqrt(var);

	// return;

	// OLD STUFF HERE
	// double mean = this->lhs->pce_coeffs[curr_timestamp][0];
	// this->pce_coeffs[curr_timestamp][0] = this->sum_coeffs * mean;
	
	// for(int i = 1; i < this->lhs->pce_coeffs[curr_timestamp].size(); i++){
	// 	this->pce_coeffs[curr_timestamp][i] = this->lhs->pce_coeffs[curr_timestamp][i] * std::sqrt(this->sum_sqr_coeffs);
	// }

	// return;

	////////////
	// double mean = this->lhs->pce_coeffs[curr_timestamp][0];
	// double var = 0.0;

	// for(int i = 0; i < this->lhs->pce_coeffs[curr_timestamp].size(); i++){
	// 	var  += lhs->pce_coeffs[curr_timestamp][i] * lhs->pce_coeffs[curr_timestamp][i] * this->bp_set_ptr->poly_sqr_expt[i];
	// }

	// this->pce_coeffs[curr_timestamp][0] = this->sum_coeffs * mean;

	// this->pce_coeffs[curr_timestamp][this->bp_set_ptr->get_var_idx(this->vars[curr_timestamp % this->vars.size()]->id)] = std::sqrt(3.0 * (this->sum_sqr_coeffs * var));

	// this->pce_coeffs[curr_timestamp][this->bp_set_ptr->get_var_idx(this->v->id)] = std::sqrt(3.0 * (this->sum_sqr_coeffs * var));

	// this->pce_coeffs[curr_timestamp][this->bp_set_ptr->get_var_idx(this->v->id)] = std::sqrt((this->sum_sqr_coeffs * var) / this->bp_set_ptr->poly_sqr_expt[this->bp_set_ptr->get_var_idx(this->v->id)]);
}