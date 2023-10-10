#include "../include/BasisPolys.h"

BasisPolySet::BasisPolySet(RandVarDist dist_type){
	this->num_vars = 0;
	this->order = 2;
	this->max_var_idx = -1;
	// this->set_size = 0;
}

BasisPolySet::~BasisPolySet(){
	for(int i = 0; i < var_arr.size(); i++){
		delete var_arr[i];
	}
	for(int i = 0; i < basis_polys.size(); i++){
		delete basis_polys[i];
	}
}

void BasisPolySet::add_variable(var* v){
	this->var_arr.push_back(v);
	this->num_vars++;

	if(v->id > this->max_var_idx){
		this->max_var_idx = v->id;
	}
}

void BasisPolySet::generate_polys(int order){
	// var* x_i;

	// for(int i = 0; i < num_vars; i++){
	// 	x_i = new var(1, -1, i+1);
	// 	this->var_arr.push_back(x_i);
	// }

	///////////////////////////////////////////////////
	std::clog << "Num vars: " << this->var_arr.size() << std::endl << std::endl;

	std::clog << "Generating basis polynomials.... ";
	this->order = order;
	std::vector<std::vector<polynomial*>> univariate_polys(this->num_vars);

	polynomial* p;
	for(int i = 0; i < this->var_arr.size(); i++){
		if(this->var_arr[i]->dist_type == GAUSSIAN){
			// std::cout << "USING HERMITE POLYS" << std::endl;
			univariate_polys[i] = std::vector<polynomial*>();

			p = new polynomial();
			p->m->coeff = 1.0;
			p->max_exp = 0;
			p->prev = nullptr;
			p->next = nullptr;

			univariate_polys[i].push_back(p);

			p = new polynomial();
			p->m->coeff = 1.0;
			p->m->arr.push_back(new term{.v = this->var_arr[i], .exp = 1});
			p->max_exp = 1;
			p->next = nullptr;
			p->prev = nullptr;
			p->var_ids_contained.push_back(this->var_arr[i]->id);

			univariate_polys[i].push_back(p);

			p = new polynomial();
			p->m->coeff = 1.0;
			p->m->arr.push_back(new term{.v = this->var_arr[i], .exp = 2});
			p->next = new polynomial();
			p->next->m->coeff = -1.0;
			p->next->prev = p;
			p->max_exp = 2;
			p->var_ids_contained.push_back(this->var_arr[i]->id);
			
			univariate_polys[i].push_back(p);
/*
			p = new polynomial();
			p->m->coeff = 1.0;
			p->m->arr.push_back(new term{.v = this->var_arr[i], .exp = 3});
			p->next = new polynomial();
			p->next->m->coeff = -3.0;
			p->next->m->arr.push_back(new term{.v = this->var_arr[i], .exp = 1});
			p->next->prev = p;
			p->max_exp = 2;
			
			univariate_polys[i].push_back(p);

			p = new polynomial();
			p->m->coeff = 1.0;
			p->m->arr.push_back(new term{.v = this->var_arr[i], .exp = 4});
			p->next = new polynomial();
			p->next->m->coeff = -6.0;
			p->next->m->arr.push_back(new term{.v = this->var_arr[i], .exp = 2});
			p->next->prev = p;
			p->next->next = new polynomial();
			p->next->next->m->coeff = 3.0;
			p->next->next->prev = p;
			p->max_exp = 2;
			
			univariate_polys[i].push_back(p);
*/
		}

		else if(this->var_arr[i]->dist_type == UNIFORM){
			// std::cout << "USING LEGENDRE POLYS" << std::endl;
			univariate_polys[i] = std::vector<polynomial*>();

			p = new polynomial();
			p->m->coeff = 1.0;
			p->max_exp = 0;
			p->prev = nullptr;
			p->next = nullptr;
			p->max_exp = 0;

			univariate_polys[i].push_back(p);

			p = new polynomial();
			p->m->coeff = 1.0;
			p->m->arr.push_back(new term{.v = this->var_arr[i], .exp = 1});
			p->max_exp = 1;
			p->next = nullptr;
			p->prev = nullptr;
			p->max_exp = 1;
			p->var_ids_contained.push_back(this->var_arr[i]->id);

			univariate_polys[i].push_back(p);

			p = new polynomial();
			p->m->coeff = 1.5;
			p->m->arr.push_back(new term{.v = this->var_arr[i], .exp = 2});
			p->next = new polynomial();
			p->next->m->coeff = -0.5;
			p->next->prev = p;
			p->max_exp = 2;
			p->var_ids_contained.push_back(this->var_arr[i]->id);
			
			univariate_polys[i].push_back(p);
	/*
			p = new polynomial();
			p->m->coeff = 2.5;
			p->m->arr.push_back(new term{.v = this->var_arr[i], .exp = 3});
			p->next = new polynomial();
			p->next->m->coeff = -1.5;
			p->next->m->arr.push_back(new term{.v = this->var_arr[i], .exp = 1});
			p->next->prev = p;
			p->max_exp = 3;
			
			univariate_polys[i].push_back(p);
	*/
		}
	}

	///////////////////////////////////////////////////

	if(this->num_vars != 0){
		faster_tensor_prod(univariate_polys, this->basis_polys, order);
	}
	else{
		p = new polynomial();
		p->m->coeff = 1.0;
		p->max_exp = 0;
		p->prev = nullptr;
		p->next = nullptr;
		p->max_exp = 0;

		this->basis_polys.push_back(p);
	}
	
	this->set_size = this->basis_polys.size();


	// TODO: Fix this very bad hack to keep accomodate else block within faster_tensor_prod
	if(univariate_polys.size() > 1){
		for(int i = 0; i < univariate_polys.size(); i++){
			for(int j = 0; j < univariate_polys[i].size(); j++){
				delete univariate_polys[i][j];
			}
		}
	}
	///////////////////////////////////////////////////
	

	gen_exp_table();
	gen_poly_expt_table();
	gen_poly_expt_sqr_table();

	std::clog << " DONE" << std::endl;
	std::clog << "    " << "Basis set has " << this->set_size << " polynomials" << std::endl;
	for(int i = 0; i < this->var_arr.size(); i++){
		if(this->var_arr[i]->dist_type == UNIFORM){
			std::cout << "    " << "ζ" << this->var_arr[i]->id << " ~ UNIFORM" << std::endl;
		}
		else if(this->var_arr[i]->dist_type == GAUSSIAN){
			std::cout << "    " << "ζ" << this->var_arr[i]->id << " ~ GAUSSIAN" << std::endl;
		}
	}

}

void BasisPolySet::regenerate_polys(int order){
	for(int i = 0; i < this->basis_polys.size(); i++){
		delete this->basis_polys[i];
	}
	this->basis_polys.clear();
	std::clog << "-----------------------" << std::endl;
	this->generate_polys(order);
}

void BasisPolySet::gen_exp_table(){
	polynomial* p1;
	polynomial* p2;
	double exp_val = 0.0;

	int n = basis_polys.size();
	this->expt_table = std::vector<std::vector<std::vector<double>>>(n, std::vector<std::vector<double>>(n, std::vector<double>(n, 0.0)));

	for(int i = 0; i < this->set_size; i++){
		for(int j = i; j < this->set_size; j++){
			for(int k = j; k < this->set_size; k++){
				p1 = new polynomial();
				p2 = new polynomial();

				mult_poly(basis_polys[i], basis_polys[j], p1);
				mult_poly(p1, basis_polys[k], p2);
				
				// exp_val = expect_poly(p2, this->dist_t);
				exp_val = expect_poly(p2);

				this->expt_table[i][j][k] = exp_val;
				this->expt_table[i][k][j] = exp_val;
				this->expt_table[j][i][k] = exp_val;
				this->expt_table[j][k][i] = exp_val;
				this->expt_table[k][i][j] = exp_val;
				this->expt_table[k][j][i] = exp_val;
				
				delete p1;
				delete p2;
			}
		}
	}
}

void BasisPolySet::gen_poly_expt_sqr_table(){
	polynomial* p;
	this->poly_sqr_expt = std::vector<double>(basis_polys.size(), 0.0);

	for(int i = 0; i < this->set_size; i++){
		p = new polynomial();
		mult_poly(this->basis_polys[i], this->basis_polys[i], p);
		// this->poly_sqr_expt[i] = expect_poly(p, this->dist_t);
		this->poly_sqr_expt[i] = expect_poly(p);
		delete p;		
	}
}

void BasisPolySet::gen_poly_expt_table(){
	this->poly_expt = std::vector<double>(basis_polys.size(), 0.0);

	for(int i = 0; i < this->set_size; i++){
		// this->poly_expt[i] = expect_poly(this->basis_polys[i], this->dist_t);
		this->poly_expt[i] = expect_poly(this->basis_polys[i]);
	}
}

int BasisPolySet::get_new_var_id(){
	if(this->max_var_idx + 1 == 0){
		this->max_var_idx = 1;
	}
	else{
		this->max_var_idx++;
	}
	
	return this->max_var_idx;
}

int BasisPolySet::get_var_idx(int id){
	for(int i = 0; i < this->basis_polys.size(); i++){
		if(this->basis_polys[i]->max_exp == 1){
			if(this->basis_polys[i]->m->arr[0]->v->id == id){
				return i;
			}
		}
	}

	return -1;
}

int BasisPolySet::get_var_idx(int id1, int id2){
	for(int i = 0; i < this->basis_polys.size(); i++){
		if(this->basis_polys[i]->var_ids_contained.size() == 2){
			if( (this->basis_polys[i]->var_ids_contained[0] == id1 || this->basis_polys[i]->var_ids_contained[0] == id2) &&
			    (this->basis_polys[i]->var_ids_contained[1] == id1 || this->basis_polys[i]->var_ids_contained[1] == id2)){
					return i;
				}
		}
	}
}

void BasisPolySet::print(){
	polynomial* p;

	for(int i = 0; i < this->set_size; i++)	{
		std::clog << "P" << i << ": ";
		this->basis_polys[i]->print();
		std::clog << std::endl;

		std::clog << "E[P" << i << "]: ";
		std::clog << expect_poly(this->basis_polys[i]);
		std::clog << std::endl;

		std::clog << "E[(P" << i << ")^2]: ";
		std::clog << poly_sqr_expt[i] << std::endl;
		
		std::clog << std::endl;
	}
}

void BasisPolySet::print_exp_table(){
	for(int i = 0; i < this->expt_table.size(); i++){
		for(int j = i; j < this->expt_table[i].size(); j++){
			for(int k = j; k < this->expt_table[j].size(); k++){
				std::cout << i << " " << j << " " << k << ": " << this->expt_table[i][j][k] << std::endl;
			}
		}
	}
}

void BasisPolySet::print_polys(){
	for(int i = 0; i < this->basis_polys.size(); i++){
		std::cout << "P" << i << ": ";
		this->basis_polys[i]->print();
		std::cout << std::endl;
	}
}