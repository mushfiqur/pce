#include "../include/dsp_block.h"

sine_node::sine_node(BasisPolySet* bp_set, std::string label) : dfg_node(SINE_BLOCK, bp_set, label) { 

}

void sine_node::init(dfg_node* arg_node){
	arg_node->tail->add_next_node(this);
	/////

	const_node* center = new const_node(this->bp_set_ptr, "center");
	
	const_node* top_const = new const_node(this->bp_set_ptr, "top-const");
	const_node* mid_const = new const_node(this->bp_set_ptr, "mid-const");
	const_node* bot_const = new const_node(this->bp_set_ptr, "bot-const");

	const_node* end_const = new const_node(this->bp_set_ptr, "end-const");

	sub_node* d = new sub_node(this->bp_set_ptr, "sin_{in}");
	mult_node* top_a = new mult_node(this->bp_set_ptr, "top_a");
	mult_node* top_b = new mult_node(this->bp_set_ptr, "top_b");
	mult_node* top_c = new mult_node(this->bp_set_ptr, "top_c");
	mult_node* mid_a = new mult_node(this->bp_set_ptr, "mid_a");
	mult_node* bot_a = new mult_node(this->bp_set_ptr, "bot_a");
	add_node* add_a = new add_node(this->bp_set_ptr, "add_a");
	add_node* add_b = new add_node(this->bp_set_ptr, "add_b");
	add_node* out = new add_node(this->bp_set_ptr, "sin(" + arg_node->label + ")");
	
	this->tail = out;

	for(int i = 0; i < this->next_nodes.size(); i++){
		this->tail->add_next_node(this->next_nodes[i]);

		if(this->next_nodes[i]->lhs == this){
			this->next_nodes[i]->lhs = this->tail;
		}

		if(this->next_nodes[i]->rhs == this){
			this->next_nodes[i]->rhs = this->tail;
		}

		for(int j = 0; j < this->next_nodes[i]->prev_nodes.size(); j++){
			if(this->next_nodes[i]->prev_nodes[j]->label == this->label){
				this->next_nodes[i]->prev_nodes.erase(this->next_nodes[i]->prev_nodes.begin() + j);
			}
		}
	}
	this->next_nodes.clear();
	
	this->head = d;
	this->head->lhs = arg_node->tail;
	this->head->rhs = nullptr;


	d->init(arg_node->tail, center);

	for(int i = 0; i < arg_node->tail->next_nodes.size(); i++){
		if(arg_node->tail->next_nodes[i]->label == d->label){
			arg_node->tail->next_nodes.erase(arg_node->tail->next_nodes.begin() + i);
			break;
		}
	}

	top_a->init(d, d);

	top_b->init(top_a, d);

	top_c->init(top_b, top_const);

	mid_a->init(top_a, mid_const);

	bot_a->init(d, bot_const);

	add_a->init(top_c, mid_a);

	add_b->init(add_a, bot_a);

	out->init(add_b, end_const);

	///// Push in the order of execution
	constants.push_back(center);
	constants.push_back(top_const);
	constants.push_back(mid_const);
	constants.push_back(bot_const);
	constants.push_back(end_const);

	nodes.push_back(d);
	nodes.push_back(top_a);
	nodes.push_back(top_b);
	nodes.push_back(top_c);
	nodes.push_back(mid_a);
	nodes.push_back(bot_a);
	nodes.push_back(add_a);
	nodes.push_back(add_b);
	nodes.push_back(out);
	

}

void sine_node::process(int curr_timestamp){
	this->last_exec_time = curr_timestamp;

	double center = this->head->lhs->pce_coeffs[curr_timestamp][0];

	this->constants[0]->pce_coeffs[curr_timestamp][0] = center;	// Set center
	this->constants[1]->pce_coeffs[curr_timestamp][0] = (1.0/6.0)*std::cos(center);		// Set top_const
	this->constants[2]->pce_coeffs[curr_timestamp][0] = (-1.0/2.0)*std::sin(center);		// Set mid_const
	this->constants[3]->pce_coeffs[curr_timestamp][0] = std::cos(center);		// Set bot_const
	this->constants[4]->pce_coeffs[curr_timestamp][0] = std::sin(center);		// Set end_const

	for(int i = 0; i < this->constants.size(); i++){
		this->constants[i]->process(curr_timestamp);
	}

	for(int i = 0; i < this->nodes.size(); i++){
		this->nodes[i]->process(curr_timestamp);
	}
}

void sine_node::set_sim_params(int tot_sim_steps, int mc_samples, int basis_set_size, SimType sim_type){
	this->last_exec_time = -1;
	this->sim_type = sim_type;
	
	for(int i = 0; i < this->constants.size(); i++){
		this->constants[i]->set_sim_params(tot_sim_steps, mc_samples, basis_set_size, sim_type);
	}

	for(int i = 0; i < this->nodes.size(); i++){
		this->nodes[i]->set_sim_params(tot_sim_steps, mc_samples, basis_set_size, sim_type);
	}
}

void sine_node::print(bool print_last){
	this->tail->print(print_last);
}

void sine_node::set_bitwidth(int width){
	// Not supported
	return;
}

void sine_node::reorder_signal_polys(){
	for(int i = 0; i < this->constants.size(); i++){
		this->constants[i]->reorder_signal_polys();
	}
	for(int i = 0; i < this->nodes.size(); i++){
		this->nodes[i]->reorder_signal_polys();
	}
}

void sine_node::save_signal_polys(){
	for(int i = 0; i < this->constants.size(); i++){
		this->constants[i]->save_signal_polys();
	}
	for(int i = 0; i < this->nodes.size(); i++){
		this->nodes[i]->save_signal_polys();
	}	

	return;
}

void sine_node::remove_signal_component(){
	for(int i = 0; i < this->constants.size(); i++){
		this->constants[i]->remove_signal_component();
	}
	for(int i = 0; i < this->nodes.size(); i++){
		this->nodes[i]->remove_signal_component();
	}	

	return;
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

sine_node::~sine_node(){
	for(int i = 0; i < this->constants.size(); i++){
		delete this->constants[i];
	}
	for(int i = 0; i < this->nodes.size(); i++){
		delete this->nodes[i];
	}
}

////////////////////////////////////////////////////////////////////////////////////

cosine_node::cosine_node(BasisPolySet* bp_set, std::string label) : dfg_node(COSINE_BLOCK, bp_set, label) { 

}

void cosine_node::init(dfg_node* arg_node){
	arg_node->tail->add_next_node(this);
	/////

	const_node* center = new const_node(this->bp_set_ptr, "center");
	
	const_node* top_const = new const_node(this->bp_set_ptr, "top-const");
	const_node* mid_const = new const_node(this->bp_set_ptr, "mid-const");
	const_node* bot_const = new const_node(this->bp_set_ptr, "bot-const");

	const_node* end_const = new const_node(this->bp_set_ptr, "end-const");

	sub_node* d = new sub_node(this->bp_set_ptr, "cos_{in}");
	mult_node* top_a = new mult_node(this->bp_set_ptr, "top_a");
	mult_node* top_b = new mult_node(this->bp_set_ptr, "top_b");
	mult_node* top_c = new mult_node(this->bp_set_ptr, "top_c");
	mult_node* mid_a = new mult_node(this->bp_set_ptr, "mid_a");
	mult_node* bot_a = new mult_node(this->bp_set_ptr, "bot_a");
	add_node* add_a = new add_node(this->bp_set_ptr, "add_a");
	add_node* add_b = new add_node(this->bp_set_ptr, "add_b");
	add_node* out = new add_node(this->bp_set_ptr, "cos(" + arg_node->label + ")");
	
	this->tail = out;

	for(int i = 0; i < this->next_nodes.size(); i++){
		this->tail->add_next_node(this->next_nodes[i]);

		if(this->next_nodes[i]->lhs == this){
			this->next_nodes[i]->lhs = this->tail;
		}

		if(this->next_nodes[i]->rhs == this){
			this->next_nodes[i]->rhs = this->tail;
		}

		for(int j = 0; j < this->next_nodes[i]->prev_nodes.size(); j++){
			if(this->next_nodes[i]->prev_nodes[j]->label == this->label){
				this->next_nodes[i]->prev_nodes.erase(this->next_nodes[i]->prev_nodes.begin() + j);
			}
		}
	}
	this->next_nodes.clear();
	
	this->head = d;
	this->head->lhs = arg_node->tail;
	this->head->rhs = nullptr;


	d->init(arg_node->tail, center);

	for(int i = 0; i < arg_node->tail->next_nodes.size(); i++){
		if(arg_node->tail->next_nodes[i]->label == d->label){
			arg_node->tail->next_nodes.erase(arg_node->tail->next_nodes.begin() + i);
			break;
		}
	}

	top_a->init(d, d);

	top_b->init(top_a, d);

	top_c->init(top_b, top_const);

	mid_a->init(top_a, mid_const);

	bot_a->init(d, bot_const);

	add_a->init(top_c, mid_a);

	add_b->init(add_a, bot_a);

	out->init(add_b, end_const);

	///// Push in the order of execution
	constants.push_back(center);
	constants.push_back(top_const);
	constants.push_back(mid_const);
	constants.push_back(bot_const);
	constants.push_back(end_const);

	nodes.push_back(d);
	nodes.push_back(top_a);
	nodes.push_back(top_b);
	nodes.push_back(top_c);
	nodes.push_back(mid_a);
	nodes.push_back(bot_a);
	nodes.push_back(add_a);
	nodes.push_back(add_b);
	nodes.push_back(out);
	

}

void cosine_node::process(int curr_timestamp){
	this->last_exec_time = curr_timestamp;

	double center = this->head->lhs->pce_coeffs[curr_timestamp][0];

	this->constants[0]->pce_coeffs[curr_timestamp][0] = center;	// Set center
	this->constants[1]->pce_coeffs[curr_timestamp][0] = (1.0/6.0)*std::sin(center);		// Set top_const
	this->constants[2]->pce_coeffs[curr_timestamp][0] = (1.0/2.0)*std::cos(center);		// Set mid_const
	this->constants[3]->pce_coeffs[curr_timestamp][0] = -1.0*std::sin(center);		// Set bot_const
	this->constants[4]->pce_coeffs[curr_timestamp][0] = std::cos(center);		// Set end_const

	for(int i = 0; i < this->constants.size(); i++){
		this->constants[i]->process(curr_timestamp);
	}

	for(int i = 0; i < this->nodes.size(); i++){
		this->nodes[i]->process(curr_timestamp);
	}
}

void cosine_node::set_sim_params(int tot_sim_steps, int mc_samples, int basis_set_size, SimType sim_type){
	this->last_exec_time = -1;
	this->sim_type = sim_type;
	
	for(int i = 0; i < this->constants.size(); i++){
		this->constants[i]->set_sim_params(tot_sim_steps, mc_samples, basis_set_size, sim_type);
	}

	for(int i = 0; i < this->nodes.size(); i++){
		this->nodes[i]->set_sim_params(tot_sim_steps, mc_samples, basis_set_size, sim_type);
	}
}

void cosine_node::print(bool print_last){
	this->tail->print(print_last);
}

void cosine_node::set_bitwidth(int width){
	// Not supported
	return;
}

void cosine_node::reorder_signal_polys(){
	for(int i = 0; i < this->constants.size(); i++){
		this->constants[i]->reorder_signal_polys();
	}
	for(int i = 0; i < this->nodes.size(); i++){
		this->nodes[i]->reorder_signal_polys();
	}
}

void cosine_node::save_signal_polys(){
	for(int i = 0; i < this->constants.size(); i++){
		this->constants[i]->save_signal_polys();
	}
	for(int i = 0; i < this->nodes.size(); i++){
		this->nodes[i]->save_signal_polys();
	}	

	return;
}

void cosine_node::remove_signal_component(){
	for(int i = 0; i < this->constants.size(); i++){
		this->constants[i]->remove_signal_component();
	}
	for(int i = 0; i < this->nodes.size(); i++){
		this->nodes[i]->remove_signal_component();
	}	

	return;
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

cosine_node::~cosine_node(){
	for(int i = 0; i < this->constants.size(); i++){
		delete this->constants[i];
	}
	for(int i = 0; i < this->nodes.size(); i++){
		delete this->nodes[i];
	}
}

////////////////////////////////////////////////////////////////////////////////////

collate_node::collate_node(BasisPolySet* bp_set, std::string label) : dfg_node(DELAY, bp_set, label){
	var* v = new var(1, -1, this->bp_set_ptr->get_new_var_id());
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
	}
	
	for(int i = 0; i < filt_coeffs.size(); i++){
		for(int j = 0; j < filt_coeffs.size(); j++){
			if(std::abs(i - j) <= correlation_dist){
				this->sum_sqr_coeffs += std::abs(filt_coeffs[i] * filt_coeffs[j]);
			}
			// else{
			// 	this->sum_sqr_coeffs += correlation * filt_coeffs[i] * filt_coeffs[j];
			// }
			
		}
	}
}

void fir_node::init(dfg_node* arg_node, std::vector<double>& coeffs){
	this->head->lhs = arg_node->tail;
	this->head->rhs = nullptr;

	arg_node->tail->add_next_node(this);

	// for(int i = 0; i < coeffs.size(); i++){
	// 	this->filt_coeffs.push_back(coeffs[i]);

	// 	this->sum_coeffs += coeffs[i];
	// 	this->sum_sqr_coeffs += (coeffs[i] * coeffs[i]);
	// }

	for(int i = 0; i < coeffs.size(); i++){

		for(int j = 0; j < coeffs.size(); j++){
			this->sum_sqr_coeffs += coeffs[i] * coeffs[j];
		}
	}
}

void fir_node::set_sim_params(int tot_sim_steps, int mc_samples, int basis_set_size, SimType sim_type){
	this->sim_type = sim_type;
	this->last_exec_time = -1;

	this->pce_coeffs = std::vector<std::vector<double>>(tot_sim_steps, std::vector<double>(basis_set_size, 0.0));

	
	if(this->tail != this){
		dfg_node* n = this->tail->rhs;	// noise node

		this->tail->pce_coeffs = std::vector<std::vector<double>>(tot_sim_steps, std::vector<double>(basis_set_size, 0.0));
		n->pce_coeffs = std::vector<std::vector<double>>(tot_sim_steps, std::vector<double>(basis_set_size, 0.0));

		for(int t = 0; t < n->pce_coeffs.size(); t++){
			n->pce_coeffs[t][this->bp_set_ptr->get_var_idx(n->node_id)] = (this->filt_coeffs.size()) * std::pow(2.0, -1.0*this->bitwidth) / 2.0;		
		}
	}
}

void fir_node::set_bitwidth(int width){
	dfg_node::set_bitwidth(width);
}

void fir_node::print(bool print_last){
	this->tail->print(print_last);
}

void fir_node::add_dist(int num_rand_vars){
	var* v;
	for(int i = 0; i < num_rand_vars; i++){
		v = new var(1, -1, this->bp_set_ptr->get_new_var_id());
		this->bp_set_ptr->add_variable(v);
		this->vars.push_back(v);
	}
	
}

void fir_node::process(int curr_timestamp){
	this->last_exec_time = curr_timestamp;

	double mean = this->lhs->pce_coeffs[curr_timestamp][0];
	this->pce_coeffs[curr_timestamp][0] = this->sum_coeffs * mean;
	
	for(int i = 1; i < this->lhs->pce_coeffs[curr_timestamp].size(); i++){
		this->pce_coeffs[curr_timestamp][i] = this->lhs->pce_coeffs[curr_timestamp][i] * std::sqrt(this->sum_sqr_coeffs);
	}

	return;

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