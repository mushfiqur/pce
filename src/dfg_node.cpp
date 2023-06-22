#include "../include/dfg_node.h"

dfg_node::dfg_node(NodeType t, std::string label, int node_id) {
	this->t = t;
	// this->pce_coeffs = std::vector<double>(set_size, 0.0);
	// this->pce_coeffs = std::vector<std::vector<double>>(tot_sim_steps, std::vector<double>(set_size, 0.0));
	// this->real_vals = std::vector<double>(tot_sim_steps, 0.0);
	this->node_id = node_id;
	this->lhs = nullptr;
	this->rhs = nullptr;
	this->last_exec_time = -1;
	this->label = label;
}

void dfg_node::process(int curr_timestamp){
	this->last_exec_time = curr_timestamp;
}

void dfg_node::init(){

}

void dfg_node::add_next_node(dfg_node * n){
	bool found = false;
	for(int i = 0; i < this->next_nodes.size(); i++){
		if(this->next_nodes[i] == n){
			found = true;
			break;
		}
	}
	if(!found){
		this->next_nodes.push_back(n);
	}

	found = false;
	for(int i = 0; i < n->prev_nodes.size(); i++){
		if(n->prev_nodes[i] == this){
			found = true;
			break;
		}
	}
	if(!found){
		n->prev_nodes.push_back(this);
	}
}

void dfg_node::print(BasisPolySet& bp_set){
	for(int t = 0; t < this->pce_coeffs.size(); t++){
		std::cerr << "[" << this->label.c_str() << " @ " << t << "] " << std::endl;
		for(int i = 0; i < this->pce_coeffs[t].size(); i++){
			if(this->pce_coeffs[t][i] != 0){
				std::cerr << this->pce_coeffs[t][i] << ": ";
				bp_set.basis_polys[i]->print();
				std::cerr << std::endl;
			}
			else{
				// std::cout << this->pce_coeffs[t][i] << std::endl;
			}
		}

		std::cerr << std::endl;
	}

}

bool dfg_node::node_args_ready(int curr_timestamp){
	bool lhs_ready = false;
	bool rhs_ready = false;
	if(this->lhs == nullptr && this->rhs == nullptr){
		return true;
	}

	if(this->lhs != nullptr) {
		if(this->lhs->last_exec_time == curr_timestamp){
			lhs_ready = true;
		}
	}
	else{
		lhs_ready = true;
	}

	if(this->rhs != nullptr) {
		if(this->rhs->last_exec_time == curr_timestamp){
			rhs_ready = true;
		}
	}
	else{
		rhs_ready = true;
	}

	return (lhs_ready && rhs_ready);

	// if((this->lhs == nullptr && this->rhs == nullptr) || 
	// ((this->lhs != nullptr && (this->lhs->last_exec_time == curr_timestamp)) && 
	// 	this->rhs != nullptr && (this->rhs->last_exec_time == curr_timestamp)))
	// {
	// 	return true;
	// }
	// else{
	// 	return false;
	// }
}

void dfg_node::print_pwr(BasisPolySet& bp_set){
	double pwr;

	for(int t = 0; t < this->pce_coeffs.size(); t++){
		pwr = 0.0;
		for(int i = 0; i < this->pce_coeffs[t].size(); i++){
			pwr += this->pce_coeffs[t][i] * this->pce_coeffs[t][i] * bp_set.poly_sqr_expt[i];
		}
		std::cerr << "[" << this->label.c_str() << " @ " << t << "] " << pwr << std::endl;
	}

}

void dfg_node::get_mc_stats(std::vector<double>& mean_arr, std::vector<double>& var_arr){
	for(int t = 0; t < this->mc_samples.size(); t++){
		mean_arr[t] = 0.0;
		for(int i = 0; i < this->mc_samples[t].size(); i++){
			mean_arr[t] += this->mc_samples[t][i];
		}
		mean_arr[t] /= this->mc_samples[t].size();
	}

	for(int t = 0; t < this->mc_samples.size(); t++){
		var_arr[t] = 0.0;
		for(int i = 0; i < this->mc_samples[t].size(); i++){
			var_arr[t] += std::pow(this->mc_samples[t][i] - mean_arr[t], 2.0);
		}

		var_arr[t] /= (this->mc_samples[t].size());
	}
}

void dfg_node::get_pce_stats(std::vector<double>& mean_arr, std::vector<double>& var_arr, BasisPolySet& bp_set){
	for(int t = 0; t < var_arr.size(); t++){
		var_arr[t] = 0.0;
		mean_arr[t] = 0.0;
		for(int i = 0; i < this->pce_coeffs[t].size(); i++){
			var_arr[t] += this->pce_coeffs[t][i] * this->pce_coeffs[t][i] * bp_set.poly_sqr_expt[i];
			mean_arr[t] += this->pce_coeffs[t][i] * bp_set.poly_expt[i];
		}
	}
}


void dfg_node::set_sim_params(int tot_sim_steps, int mc_samples, int basis_set_size, SimType sim_type){
	this->sim_type = sim_type;
	this->last_exec_time = -1;

	if(sim_type == PCE){
		this->pce_coeffs = std::vector<std::vector<double>>(tot_sim_steps, std::vector<double>(basis_set_size, 0.0));
	}
	else{
		this->mc_samples = std::vector<std::vector<double>>(tot_sim_steps, std::vector<double>(mc_samples, 0.0));
	}

	if(this->t == INPUT){
		input_node* n;
		n = static_cast<input_node*>(this);
		n->set_sim_params();
	}
	else if(this->t == CONST){
		const_node* n;
		n = static_cast<const_node*>(this);
		n->set_sim_params();
	}
}

///////////////////////////////////////////////////

input_node::input_node(BasisPolySet& bp_set, int node_id, std::string label) : dfg_node(INPUT, label, node_id){
	bool var_exists = false;
	for(int i = 0; i < bp_set.var_arr.size(); i++){
		if(bp_set.var_arr[i]->id == node_id){
			var_exists = true;
			this->v = bp_set.var_arr[i];
			break;
		}
	}

	if(!var_exists){
		var* v = new var(1, -1, node_id);
		bp_set.add_variable(v);
		this->v = v;
	}
	
	std::random_device rd;
	this->mt = std::mt19937(rd());
	this->dist = std::uniform_real_distribution<double>(this->v->a, this->v->b);
};


void input_node::init(){
	// this->sim_type = sim_type;

	// for(int i = 0; i < this->pce_coeffs.size(); i++){
	// 	this->pce_coeffs[i][idx] = 1.0;

	// }
}

void input_node::set_sim_params(){
	if(sim_type == PCE){
		for(int t = 0; t < this->pce_coeffs.size(); t++){
			this->pce_coeffs[t][this->v->id] = 1.0;
		}
	}
	else{
		std::vector<double> samples(this->mc_samples[0].size(), 0.0);

		for(int i = 0; i < samples.size(); i++)		{
			samples[i] = this->dist(this->mt);
		}

		for(int i = 0; i < this->mc_samples.size(); i++){
			this->mc_samples[i] = samples;
		}
	}
}


void input_node::process(int curr_timestamp){
	this->last_exec_time = curr_timestamp;

	// if(this->sim_type == MONTE_CARLO){
	// 	this->real_vals[curr_timestamp] = this->dist(this->mt);
	// }
}

///////////////////////////////////////////////////

mult_node::mult_node(BasisPolySet& bp_set, std::string label) : dfg_node(MULT, label){
	this->ref_exp_table = &bp_set.expt_table;
	this->ref_exp_sqr_table = &bp_set.poly_sqr_expt;
};

void mult_node::init(dfg_node* lhs, dfg_node* rhs){
	this->lhs = lhs;
	this->rhs = rhs;

	lhs->add_next_node(this);
	rhs->add_next_node(this);
}

void mult_node::process(int curr_timestamp){
	this->last_exec_time = curr_timestamp;

	if(this->sim_type == PCE){
		process_pce_sim(curr_timestamp);
	}else{
		process_mc_sim(curr_timestamp);
	}
}

void mult_node::process_pce_sim(int curr_timestamp){
	if(lhs->t == CONST || rhs->t == CONST){
		// Linear Multiplication
		const_node* c;
		dfg_node* n;
		if(lhs->t == CONST){
			c = static_cast<const_node*>(lhs);
			n = rhs;
		}
		else{
			c = static_cast<const_node*>(rhs);
			n = lhs;
		}

		for(int i = 0; i < n->pce_coeffs[curr_timestamp].size(); i++){
			this->pce_coeffs[curr_timestamp][i] = c->pce_coeffs[curr_timestamp][0] * n->pce_coeffs[curr_timestamp][i];
		}
	}
	else{
		// Non-linear Multiplication
		for(int k = 0; k < this->pce_coeffs[curr_timestamp].size(); k++){
			for(int i = 0; i < lhs->pce_coeffs[curr_timestamp].size(); i++){
				for(int j = 0; j < rhs->pce_coeffs[curr_timestamp].size(); j++){
					this->pce_coeffs[curr_timestamp][k] += lhs->pce_coeffs[curr_timestamp][i] * rhs->pce_coeffs[curr_timestamp][j] * (*ref_exp_table)[i][j][k];
				}
			}
			this->pce_coeffs[curr_timestamp][k] /= (*ref_exp_sqr_table)[k];
		}
	}
}

void mult_node::process_mc_sim(int curr_timestamp){
	if(lhs->t == CONST || rhs->t == CONST){
		const_node* c;
		dfg_node* n;
		if(lhs->t == CONST){
			c = static_cast<const_node*>(lhs);
			n = rhs;
		}
		else{
			c = static_cast<const_node*>(rhs);
			n = lhs;
		}

		for(int i = 0; i < n->mc_samples[curr_timestamp].size(); i++){
			this->mc_samples[curr_timestamp][i] = c->mc_samples[curr_timestamp][0] * n->mc_samples[curr_timestamp][i];
		}

	}
	else{
		for(int i = 0; i < this->mc_samples[curr_timestamp].size(); i++){
			this->mc_samples[curr_timestamp][i] = lhs->mc_samples[curr_timestamp][i] * rhs->mc_samples[curr_timestamp][i];
		}
	}
	// this->real_vals[curr_timestamp] = lhs->real_vals[curr_timestamp] * rhs->real_vals[curr_timestamp];
}

///////////////////////////////////////////////////

void add_node::init(dfg_node* lhs, dfg_node* rhs){
	this->lhs = lhs;
	this->rhs = rhs;

	lhs->add_next_node(this);
	rhs->add_next_node(this);
}

void add_node::process(int curr_timestamp){
	this->last_exec_time = curr_timestamp;

	if(this->sim_type == PCE){
		process_pce_sim(curr_timestamp);
	}
	else{
		process_mc_sim(curr_timestamp);
	}
}

void add_node::process_mc_sim(int curr_timestamp){
	for(int i = 0; i < this->mc_samples[curr_timestamp].size(); i++){
		this->mc_samples[curr_timestamp][i] = lhs->mc_samples[curr_timestamp][i] + rhs->mc_samples[curr_timestamp][i];
	}
}

void add_node::process_pce_sim(int curr_timestamp){
	for(int i = 0; i < this->pce_coeffs[curr_timestamp].size(); i++){
		this->pce_coeffs[curr_timestamp][i] = lhs->pce_coeffs[curr_timestamp][i] + rhs->pce_coeffs[curr_timestamp][i];
	}
}

///////////////////////////////////////////////////

void sub_node::init(dfg_node* lhs, dfg_node* rhs){
	this->lhs = lhs;
	this->rhs = rhs;

	lhs->add_next_node(this);
	rhs->add_next_node(this);
}

void sub_node::process(int curr_timestamp){
	this->last_exec_time = curr_timestamp;
	
	if(this->sim_type == PCE){
		process_pce_sim(curr_timestamp);
	}
	else{
		process_mc_sim(curr_timestamp);
	}

}

void sub_node::process_pce_sim(int curr_timestamp){
	for(int i = 0; i < this->pce_coeffs.size(); i++){
		this->pce_coeffs[curr_timestamp][i] = lhs->pce_coeffs[curr_timestamp][i] - rhs->pce_coeffs[curr_timestamp][i];
	}
}

void sub_node::process_mc_sim(int curr_timestamp){
	for(int i = 0; i < this->mc_samples[curr_timestamp].size(); i++){
		this->mc_samples[curr_timestamp][i] = lhs->mc_samples[curr_timestamp][i] - rhs->mc_samples[curr_timestamp][i];
	}
}

///////////////////////////////////////////////////

void delay_node::init(dfg_node* prev_node){
	this->lhs = prev_node;
	this->rhs = nullptr;

	this->lhs->add_next_node(this);
}

void delay_node::process(int curr_timestamp){
	this->last_exec_time = curr_timestamp;

	if(this->sim_type == PCE){
		process_pce_sim(curr_timestamp);
	}
	else{
		process_mc_sim(curr_timestamp);
	}
}

void delay_node::process_pce_sim(int curr_timestamp){
	if(curr_timestamp - 1 < 0){
		for(int i = 0; i < this->pce_coeffs[curr_timestamp].size(); i++){
			this->pce_coeffs[curr_timestamp][i] = 0.0;
		}
	}
	else{
		for(int i = 0; i < lhs->pce_coeffs[curr_timestamp-1].size(); i++){
			this->pce_coeffs[curr_timestamp][i] = lhs->pce_coeffs[curr_timestamp - 1][i];
		}
	}
}

void delay_node::process_mc_sim(int curr_timestamp){
	if(curr_timestamp - 1 < 0){
		for(int i = 0; i < this->mc_samples[curr_timestamp].size(); i++){
			this->mc_samples[curr_timestamp][i] = 0.0;
		}
	}
	else{
		for(int i = 0; i < this->mc_samples[curr_timestamp].size(); i++){
			this->mc_samples[curr_timestamp][i] = lhs->mc_samples[curr_timestamp - 1][i];
		}
	}
}

///////////////////////////////////////////////////

void const_node::init(double val){
	this->val = val;
}

void const_node::set_sim_params(){
	if(this->sim_type == PCE){
		for(int i = 0; i < pce_coeffs.size(); i++){
			this->pce_coeffs[i][0] = this->val;
		}
	}
	else{
		for(int i = 0; i < mc_samples.size(); i++){
			this->mc_samples[i][0] = this->val;
		}
	}
}

void const_node::process(int curr_timestamp){
	this->last_exec_time = curr_timestamp;
}