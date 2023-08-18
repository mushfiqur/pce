#include "../include/dfg_node.h"

dfg_node::dfg_node(NodeType t, BasisPolySet* bp_set, std::string label, int node_id) {
	this->t = t;
	this->bp_set_ptr = bp_set;
	this->node_id = node_id;
	this->lhs = nullptr;
	this->rhs = nullptr;
	this->last_exec_time = -1;
	this->label = label;
	this->bitwidth = -1;

	this->head = this;
	this->tail = this;
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

void dfg_node::print(bool print_last){
	if(print_last){
		int last_idx = this->pce_coeffs.size() - 1;

		std::cerr << "[" << this->label.c_str() << " @ " << last_idx << "] " << std::endl;
		for(int i = 0; i < this->pce_coeffs[t].size(); i++){
			if(this->pce_coeffs[last_idx][i] != 0){
				std::cerr << this->pce_coeffs[last_idx][i] << ": ";
				this->bp_set_ptr->basis_polys[i]->print();
				std::cerr << std::endl;
			}
			else{
				// std::cout << this->pce_coeffs[t][i] << std::endl;
			}
		}

		std::cerr << std::endl;
	}
	else{
		for(int t = 0; t < this->pce_coeffs.size(); t++){
			std::cerr << "[" << this->label.c_str() << " @ " << t << "] " << std::endl;
			for(int i = 0; i < this->pce_coeffs[t].size(); i++){
				if(this->pce_coeffs[t][i] != 0){
					std::cerr << this->pce_coeffs[t][i] << ": ";
					this->bp_set_ptr->basis_polys[i]->print();
					std::cerr << std::endl;
				}
				else{
					// std::cout << this->pce_coeffs[t][i] << std::endl;
				}
			}

			std::cerr << std::endl;
		}
	}

}

void dfg_node::print_signal_coeffs(){
	for(int t = 0; t < this->signal_coeffs.size(); t++){
		std::cerr << "[" << this->label.c_str() << " @ " << t << "] " << std::endl;
		for(int i = 0; i < this->signal_coeffs[t].size(); i++){
			if(this->signal_coeffs[t][i] != 0){
				std::cerr << this->signal_coeffs[t][i] << ": ";
				this->bp_set_ptr->basis_polys[i]->print();
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

}

void dfg_node::set_bitwidth(int width){
	// this->bitwidth = width;

	if(this->tail == this){
		int n_noise_id = this->bp_set_ptr->get_new_var_id();
		std::string lbl = "n" + std::to_string(n_noise_id);
		noise_node* n_noise = new noise_node(this->bp_set_ptr, n_noise_id, lbl);
		n_noise->bitwidth = width;

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

void dfg_node::remove_signal_component(){
	for(int t = 0; t < this->pce_coeffs.size(); t++){
		for(int i = 0; i < this->pce_coeffs[t].size(); i++){
			this->pce_coeffs[t][i] -= this->signal_coeffs[t][i];
			if(this->tail != this){
				this->tail->pce_coeffs[t][i] -= this->signal_coeffs[t][i];
			}
		}
	}
}

double dfg_node::get_pwr(){
	double pwr = 0.0;
	int last_idx = this->pce_coeffs.size() - 1;

	for(int i = 0; i < this->pce_coeffs[last_idx].size(); i++){
		pwr += this->pce_coeffs[last_idx][i] * this->pce_coeffs[last_idx][i] * this->bp_set_ptr->poly_sqr_expt[i];
	}

	return pwr;
}

void dfg_node::print_pwr(){
	double pwr;

	for(int t = 0; t < this->pce_coeffs.size(); t++){
		pwr = 0.0;
		for(int i = 0; i < this->pce_coeffs[t].size(); i++){
			pwr += this->pce_coeffs[t][i] * this->pce_coeffs[t][i] * this->bp_set_ptr->poly_sqr_expt[i];
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
			var_arr[t] += std::pow(this->mc_samples[t][i], 2.0);
			// var_arr[t] += std::pow(this->mc_samples[t][i] - mean_arr[t], 2.0);
		}

		var_arr[t] /= (this->mc_samples[t].size());
	}
}

void dfg_node::get_pce_stats(std::vector<double>& mean_arr, std::vector<double>& var_arr){
	for(int t = 0; t < var_arr.size(); t++){
		var_arr[t] = 0.0;
		mean_arr[t] = 0.0;
		for(int i = 0; i < this->pce_coeffs[t].size(); i++){
			var_arr[t] += this->pce_coeffs[t][i] * this->pce_coeffs[t][i] * this->bp_set_ptr->poly_sqr_expt[i];
			mean_arr[t] += this->pce_coeffs[t][i] * this->bp_set_ptr->poly_expt[i];
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

	if(this->tail != this){
		this->tail->set_sim_params(tot_sim_steps, mc_samples, basis_set_size, sim_type);
		this->tail->rhs->set_sim_params(tot_sim_steps, mc_samples, basis_set_size, sim_type);
	}
}

void dfg_node::save_signal_polys(){
	this->signal = std::vector<std::vector<signal_polys>>(this->pce_coeffs.size(), std::vector<signal_polys>(this->bp_set_ptr->basis_polys.size(), {.coeff = 0.0, .poly=nullptr}));

	for(int t = 0; t < this->pce_coeffs.size(); t++){
		for(int i = 0; i < this->pce_coeffs[t].size(); i++){
			this->signal[t][i].coeff = this->pce_coeffs[t][i];
			this->signal[t][i].poly = this->bp_set_ptr->basis_polys[i]->copy();
		}
	}
}

void dfg_node::reorder_signal_polys(){
	this->signal_coeffs = std::vector<std::vector<double>>(this->pce_coeffs.size(), std::vector<double>(this->bp_set_ptr->basis_polys.size(), 0.0));
	
	// Time is but a number
	for(int t = 0; t < this->signal.size(); t++){
		for(int i = 0; i < this->signal[t].size(); i++){
			for(int j = 0; j < this->bp_set_ptr->basis_polys.size(); j++){
				if(this->signal[t][i].poly->equals(this->bp_set_ptr->basis_polys[j])){
					this->signal_coeffs[t][j] = this->signal[t][i].coeff;
				}
			}
		}
	}
}

dfg_node::~dfg_node(){
	if(this->tail != this && this->t != SINE_BLOCK && this->t != COSINE_BLOCK){
		delete this->tail->rhs;
		delete this->tail;
	}
	
	for(int t = 0; t < this->signal.size(); t++){
		for(int i = 0; i < this->signal[t].size(); i++){
			if(this->signal[t][i].poly != nullptr){
				delete this->signal[t][i].poly;
			}
		}
	}
}

///////////////////////////////////////////////////

input_node::input_node(BasisPolySet* bp_set, std::string label) : dfg_node(INPUT_SIGNAL, bp_set, label){
	this->bp_set_ptr = bp_set;
	this->set_range(-1.0, 1.0);
	
	// TODO: FIX - re-enable support for Monte Carlo sims
	// std::random_device rd;
	// this->mt = std::mt19937(rd());
	// this->dist = std::uniform_real_distribution<double>(this->v->a, this->v->b);
};


void input_node::init(){

}

void input_node::set_range(double a, double b){
	this->unfm_dist_param_a = a;
	this->unfm_dist_param_b = b;
}

void input_node::set_sim_params(int tot_sim_steps, int mc_samples, int basis_set_size, SimType sim_type){
	this->sim_type = sim_type;
	this->last_exec_time = -1;

	if(sim_type == PCE){
		this->pce_coeffs = std::vector<std::vector<double>>(tot_sim_steps, std::vector<double>(basis_set_size, 0.0));
		
		for(int t = 0; t < this->pce_coeffs.size(); t++){
			if(this->cfg.has_signal){
				if(this->cfg.waveform == COSINE){
					this->pce_coeffs[t][0] += std::cos(2.0*M_PI*(this->cfg.freq)*t);
				}
				else if(this->cfg.waveform == SINE){
					this->pce_coeffs[t][0] += std::sin(2.0*M_PI*(this->cfg.freq)*t);
				}
				else if(this->cfg.waveform == DC){
					this->pce_coeffs[t][0] += this->cfg.freq;
				}
			}

			if(this->cfg.has_dist){
				this->pce_coeffs[t][0] += (unfm_dist_param_a + unfm_dist_param_b) / 2.0;
				if(t < this->vars.size()){
					this->pce_coeffs[t][this->bp_set_ptr->get_var_idx(this->vars[t]->id)] = (unfm_dist_param_b - unfm_dist_param_a) / 2.0;
				}
				else{
					this->pce_coeffs[t][this->bp_set_ptr->get_var_idx(this->vars[ this->vars.size() - 1 ]->id)] = (unfm_dist_param_b - unfm_dist_param_a) / 2.0;
				}
			}


			// this->pce_coeffs[t][0] = (unfm_dist_param_a + unfm_dist_param_b) / 2.0;

			// this->pce_coeffs[t][this->bp_set_ptr->get_var_idx(this->vars[0]->id)] = (unfm_dist_param_b - unfm_dist_param_a) / 2.0;

			// this->pce_coeffs[t][this->bp_set_ptr->get_var_idx(this->v->id)] = (unfm_dist_param_b - unfm_dist_param_a) / 2.0;
		}
	}
	else{
		this->mc_samples = std::vector<std::vector<double>>(tot_sim_steps, std::vector<double>(mc_samples, 0.0));
		std::vector<double> samples(this->mc_samples[0].size(), 0.0);

		for(int i = 0; i < samples.size(); i++){
			samples[i] = this->dist(this->mt);
		}

		for(int i = 0; i < this->mc_samples.size(); i++){
			this->mc_samples[i] = samples;
		}
	}

	// TODO: Verify this hack
	if(this->tail != this){
		this->tail->set_sim_params(tot_sim_steps, mc_samples, basis_set_size, sim_type);
		this->tail->rhs->set_sim_params(tot_sim_steps, mc_samples, basis_set_size, sim_type);
	}
}

void input_node::process(int curr_timestamp){
	this->last_exec_time = curr_timestamp;
}

void input_node::set_bitwidth(int width){
	// this->bitwidth = width;
	dfg_node::set_bitwidth(width);
	this->dist = std::uniform_real_distribution<double>(0.0, std::pow(2.0, -1.0*width));
}

void input_node::add_signal(WaveType wave, double freq){
	this->cfg.waveform = wave;
	this->cfg.freq = freq;
	this->cfg.has_signal = true;
}

void input_node::add_dist(int num_rand_vars){
	var* v;
	for(int i = 0; i < num_rand_vars; i++){
		v = new var(1, -1, this->bp_set_ptr->get_new_var_id());
		this->bp_set_ptr->add_variable(v);
		this->vars.push_back(v);
	}

	this->cfg.has_dist = true;
	this->cfg.propagate_unique_dists = true;
}

///////////////////////////////////////////////////

noise_node::noise_node(BasisPolySet* bp_set, int node_id, std::string label) : dfg_node(INPUT_NOISE, bp_set, label, node_id){
	this->bp_set_ptr = bp_set;
	
	bool var_exists = false;
	for(int i = 0; i < bp_set->var_arr.size(); i++){
		if(bp_set->var_arr[i]->id == node_id){
			var_exists = true;
			this->v = bp_set->var_arr[i];
			break;
		}
	}

	if(!var_exists){
		var* v = new var(1, -1, node_id);
		bp_set->add_variable(v);
		this->v = v;
	}
	
	std::random_device rd;
	this->mt = std::mt19937(rd());
	this->dist = std::uniform_real_distribution<double>(this->v->a, this->v->b);
};


void noise_node::init(){

}

void noise_node::set_sim_params(int tot_sim_steps, int mc_samples, int basis_set_size, SimType sim_type){
	this->sim_type = sim_type;
	this->last_exec_time = -1;

	if(sim_type == PCE){
		this->pce_coeffs = std::vector<std::vector<double>>(tot_sim_steps, std::vector<double>(basis_set_size, 0.0));
		
		for(int t = 0; t < this->pce_coeffs.size(); t++){
			// TODO: line below uncommented represents a truncation quantizer
			//         commenting it out makes it a rounding quantizer
			//        Implement enum to differentiate
			// this->pce_coeffs[t][0] = std::pow(2.0, -1.0*this->bitwidth) / 2.0;
			this->pce_coeffs[t][this->bp_set_ptr->get_var_idx(this->v->id)] = std::pow(2.0, -1.0*this->bitwidth) / 2.0;
		}
	}
	else{
		this->mc_samples = std::vector<std::vector<double>>(tot_sim_steps, std::vector<double>(mc_samples, 0.0));
		std::vector<double> samples(this->mc_samples[0].size(), 0.0);

		for(int i = 0; i < samples.size(); i++){
			samples[i] = this->dist(this->mt);
		}

		for(int i = 0; i < this->mc_samples.size(); i++){
			this->mc_samples[i] = samples;
		}
	}

	// TODO: Verify this hack
	if(this->tail != this){
		this->tail->set_sim_params(tot_sim_steps, mc_samples, basis_set_size, sim_type);
		this->tail->rhs->set_sim_params(tot_sim_steps, mc_samples, basis_set_size, sim_type);
	}
}

void noise_node::process(int curr_timestamp){
	this->last_exec_time = curr_timestamp;
}

void noise_node::set_bitwidth(int width){
	this->bitwidth = width;
	this->dist = std::uniform_real_distribution<double>(0.0, std::pow(2.0, -1.0*width));
}

///////////////////////////////////////////////////

mult_node::mult_node(BasisPolySet* bp_set, std::string label) : dfg_node(MULT, bp_set, label){
	this->bp_set_ptr = bp_set;
};

void mult_node::init(dfg_node* lhs, dfg_node* rhs){
	this->head->lhs = lhs->tail;
	this->head->rhs = rhs->tail;

	lhs->tail->add_next_node(this->head);
	rhs->tail->add_next_node(this->head);
}

void mult_node::set_bitwidth(int width){
	dfg_node::set_bitwidth(width);
}

void mult_node::set_sim_params(int tot_sim_steps, int mc_samples, int basis_set_size, SimType sim_type){
	this->sim_type = sim_type;
	this->last_exec_time = -1;

	if(sim_type == PCE){
		this->pce_coeffs = std::vector<std::vector<double>>(tot_sim_steps, std::vector<double>(basis_set_size, 0.0));
	}
	else{
		this->mc_samples = std::vector<std::vector<double>>(tot_sim_steps, std::vector<double>(mc_samples, 0.0));
	}
	
	if(this->tail != this){
		this->tail->set_sim_params(tot_sim_steps, mc_samples, basis_set_size, sim_type);
		this->tail->rhs->set_sim_params(tot_sim_steps, mc_samples, basis_set_size, sim_type);
	}
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
					this->pce_coeffs[curr_timestamp][k] += lhs->pce_coeffs[curr_timestamp][i] * rhs->pce_coeffs[curr_timestamp][j] * this->bp_set_ptr->expt_table[i][j][k];
				}
			}
			this->pce_coeffs[curr_timestamp][k] /= this->bp_set_ptr->poly_sqr_expt[k];
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

// mult_node::~mult_node(){
// 	if(this->tail != this){
// 		delete this->tail->rhs;
// 		delete this->tail;
// 	}
// }

///////////////////////////////////////////////////

divide_node::divide_node(BasisPolySet* bp_set, std::string label) : dfg_node(DIVIDE, bp_set, label){
	this->bp_set_ptr = bp_set;
};

void divide_node::init(dfg_node* lhs, dfg_node* rhs){
	this->head->lhs = lhs->tail;
	this->head->rhs = rhs->tail;

	lhs->tail->add_next_node(this->head);
	rhs->tail->add_next_node(this->head);
}

void divide_node::process(int curr_timestamp){
	this->last_exec_time = curr_timestamp;

	if(this->sim_type == PCE){
		this->process_pce_sim(curr_timestamp);
	}else{
		this->process_mc_sim(curr_timestamp);
	}
}

void divide_node::process_pce_sim(int curr_timestamp){
	Eigen::VectorXd x;
	double limit = 1e-15;

	this->b = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(lhs->pce_coeffs[curr_timestamp].data(), lhs->pce_coeffs[curr_timestamp].size());

	for(int x_idx = 0; x_idx < this->bp_set_ptr->basis_polys.size(); x_idx++){
		for(int z_idx = 0; z_idx < this->bp_set_ptr->basis_polys.size(); z_idx++){
			A(x_idx,z_idx) = 0.0;

			for(int i = 0; i < this->bp_set_ptr->basis_polys.size(); i++){
				A(x_idx,z_idx) += (this->bp_set_ptr->expt_table[x_idx][z_idx][i] / this->bp_set_ptr->poly_sqr_expt[x_idx]) * rhs->pce_coeffs[curr_timestamp][i];
			}
		}
	}

	this->solver.compute(A);
	x = solver.solve(this->b);

	this->pce_coeffs[curr_timestamp] = std::vector<double>(x.data(), x.data() + x.rows() * x.cols());

	for(int i = 0; i < this->pce_coeffs[curr_timestamp].size(); i++){
		if(this->pce_coeffs[curr_timestamp][i] < limit){
			this->pce_coeffs[curr_timestamp][i] = 0.0;
		}
	}

}

void divide_node::process_mc_sim(int curr_timestamp){
	// TODO
	return;
}

void divide_node::set_bitwidth(int width){
	dfg_node::set_bitwidth(width);
}

void divide_node::set_sim_params(int tot_sim_steps, int mc_samples, int basis_set_size, SimType sim_type){
	this->sim_type = sim_type;
	this->last_exec_time = -1;

	if(sim_type == PCE){
		this->A.resize(basis_set_size, basis_set_size);
		this->pce_coeffs = std::vector<std::vector<double>>(tot_sim_steps, std::vector<double>(basis_set_size, 0.0));
	}
	else{
		this->mc_samples = std::vector<std::vector<double>>(tot_sim_steps, std::vector<double>(mc_samples, 0.0));
	}
	
	if(this->tail != this){
		this->tail->set_sim_params(tot_sim_steps, mc_samples, basis_set_size, sim_type);
		this->tail->rhs->set_sim_params(tot_sim_steps, mc_samples, basis_set_size, sim_type);
	}
}

///////////////////////////////////////////////////

void add_node::init(dfg_node* lhs, dfg_node* rhs){
	this->head->lhs = lhs->tail;
	this->head->rhs = rhs->tail;

	lhs->tail->add_next_node(this->head);
	rhs->tail->add_next_node(this->head);
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

void add_node::set_bitwidth(int width){
	dfg_node::set_bitwidth(width);
	// if(width == -1){
	// 	this->bitwidth = lhs->tail->bitwidth > rhs->tail->bitwidth ? lhs->tail->bitwidth : rhs->tail->bitwidth;
	// }
	// else{
	// 	this->bitwidth = width;
	// }
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
	this->head->lhs = lhs->tail;
	this->head->rhs = rhs->tail;

	lhs->tail->add_next_node(this->head);
	rhs->tail->add_next_node(this->head);
}

void sub_node::set_bitwidth(int width){
	dfg_node::set_bitwidth(width);
	// if(width == -1){
	// 	this->bitwidth = lhs->tail->bitwidth > rhs->tail->bitwidth ? lhs->tail->bitwidth : rhs->tail->bitwidth;
	// }
	// else{
	// 	this->bitwidth = width;
	// }
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
	for(int i = 0; i < this->pce_coeffs[curr_timestamp].size(); i++){
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
	this->head->lhs = prev_node->tail;
	this->head->rhs = nullptr;

	prev_node->tail->add_next_node(this);
}

void delay_node::set_bitwidth(int width){
	dfg_node::set_bitwidth(width);
	// if(width == -1){
	// 	this->bitwidth = lhs->tail->bitwidth;
	// }
	// else{
	// 	this->bitwidth = width;
	// }
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

void const_node::set_sim_params(int tot_sim_steps, int mc_samples, int basis_set_size, SimType sim_type){
	this->sim_type = sim_type;
	this->last_exec_time = -1;
	
	if(this->sim_type == PCE){
		this->pce_coeffs = std::vector<std::vector<double>>(tot_sim_steps, std::vector<double>(basis_set_size, 0.0));

		for(int i = 0; i < pce_coeffs.size(); i++){
			this->pce_coeffs[i][0] = this->val;
		}
	}
	else{
		this->mc_samples = std::vector<std::vector<double>>(tot_sim_steps, std::vector<double>(mc_samples, 0.0));
		for(int i = 0; i < this->mc_samples.size(); i++){
			this->mc_samples[i][0] = this->val;
		}
	}
}

void const_node::process(int curr_timestamp){
	this->last_exec_time = curr_timestamp;
}