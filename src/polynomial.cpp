#include "../include/BasisPolys.h"
#include "../include/utils.h"

var::var(int id){
	this->id = id;
}

void var::init_uniform_var(double a, double b){
	this->a = a <= b ? a : b;
	this->b = b > a ? b : a;
	this->dist_type = UNIFORM;
}

void var::init_gaussian_var(double mean, double variance){
	this->mean = mean;
	this->variance = variance;
	this->dist_type = GAUSSIAN;
}

polynomial::polynomial(){
	this->m = new monomial();
	this->next = nullptr;
	this->prev = nullptr;
}

polynomial::~polynomial(){
	delete this->m;

	if(this->next != nullptr){
		delete this->next;
	}
}

void polynomial::print(){
	this->m->print();
	if(this->next != nullptr){
		std::clog << " + ";
		this->next->print();
	}
}

polynomial* polynomial::copy(){
	polynomial* ret_poly = new polynomial();
	polynomial* p_itr = this;
	polynomial* r_itr = ret_poly;

	while(p_itr != nullptr){
		r_itr->m->coeff = p_itr->m->coeff;
		for(int i = 0; i < p_itr->m->arr.size(); i++){
			r_itr->m->arr.push_back(new term{ .v = p_itr->m->arr[i]->v, .exp = p_itr->m->arr[i]->exp});
		}

		if(p_itr->next != nullptr){
			r_itr->next = new polynomial();
			r_itr->next->prev = r_itr;
			r_itr = r_itr->next;
		}

		p_itr = p_itr->next;
	}

	return ret_poly;
}

bool polynomial::equals(polynomial* p){
	polynomial* curr_poly = this;
	polynomial* comp_poly = p;

	while(curr_poly != nullptr){
		if(comp_poly == nullptr || curr_poly->m->arr.size() != comp_poly->m->arr.size() || curr_poly->m->coeff != comp_poly->m->coeff){
			return false;
		}


		for(int i = 0; i < curr_poly->m->arr.size(); i++){
			if(curr_poly->m->arr[i]->v != comp_poly->m->arr[i]->v || curr_poly->m->arr[i]->exp != comp_poly->m->arr[i]->exp){
				return false;
			}
		}

		curr_poly = curr_poly->next;
		comp_poly = comp_poly->next;
	}

	return true;
}

monomial::monomial(){
	this->arr = std::vector<term*>();
}

monomial::~monomial(){
	for(int i = 0; i < this->arr.size(); i++){
		delete this->arr[i];
	}
}

void monomial::print(){
	std::clog << this->coeff;

	for(int i = 0; i < this->arr.size(); i++){
		std::clog << "(ζ" << this->arr[i]->v->id << ")^(" << this->arr[i]->exp << ")";
	}
}

void tensor_prod(std::vector<std::vector<polynomial*>>& univariate_set, std::vector<polynomial*>& result, int max_order){
	polynomial* p;

	for(int a = 0; a < univariate_set.size() - 1; a++){
		for(int b = a+1; b < univariate_set.size(); b++){
			for(int i = 0; i < univariate_set[a].size(); i++){
				for(int j = 0; j < univariate_set[b].size(); j++){
					if((univariate_set[a][i]->max_exp + univariate_set[b][j]->max_exp) <= max_order){
						p = new polynomial();
						mult_poly(univariate_set[a][i], univariate_set[b][j], p);
						result.push_back(p);
					}
				}
			}
		}
	}

	
}

void faster_tensor_prod(std::vector<std::vector<polynomial*>>& univariate_set, std::vector<polynomial*>& result, int max_order){
	if(univariate_set.size() > 1){
		std::vector<polynomial*> set_temp = std::vector<polynomial*>();

		tensor_prod(univariate_set[0], univariate_set[1], result, max_order);

		// call tensor prod with result set and next set, store res in result
		
		// Can improve memory usage by keeping storing result in result array instead of set_temp
		//	array. Keep track of end using pointer
		
		for(int i = 2; i < univariate_set.size(); i++){
			tensor_prod(result, univariate_set[i], set_temp, max_order);

			// result.clear();
			for(int j = 0; j < set_temp.size(); j++){
				if(j < result.size()){
					delete result[j];
					result[j] = set_temp[j];
				}
				else{
					result.push_back(set_temp[j]);
				}
			}
			set_temp = std::vector<polynomial*>();
		}
	}
	else{
		for(int i = 0; i < univariate_set[0].size(); i++){
			result.push_back(univariate_set[0][i]);
		}
	}

	
}


void tensor_prod(std::vector<polynomial*>& set_a, std::vector<polynomial*>& set_b, std::vector<polynomial*>& result, int max_order){
	polynomial* p;

	for(int i = 0; i < set_a.size(); i++){
		for(int j = 0; j < set_b.size(); j++){
			if((set_a[i]->max_exp + set_b[j]->max_exp) <= max_order){
				p = new polynomial();
				p->max_exp = set_a[i]->max_exp + set_b[j]->max_exp;
				mult_poly(set_a[i], set_b[j], p);
				
				for(int n = 0; n < set_a[i]->var_ids_contained.size(); n++){
					p->var_ids_contained.push_back(set_a[i]->var_ids_contained[n]);
				}
				
				for(int n = 0; n < set_b[j]->var_ids_contained.size(); n++){
					p->var_ids_contained.push_back(set_b[j]->var_ids_contained[n]);
				}
				
				result.push_back(p);
			}
		}
	}
		
	
}

void mult_poly(polynomial* p, polynomial*q, polynomial* result){
	polynomial* p_itr = p;
	polynomial* q_itr = q;
	polynomial* r_itr = result;

	// Multiply monomials
	while(p_itr != nullptr){
		while(q_itr != nullptr){
			mult_mono(p_itr->m, q_itr->m, r_itr->m);
			
			q_itr = q_itr->next;
			if(q_itr != nullptr){
				r_itr->next = new polynomial();
				r_itr->next->prev = r_itr;			
				r_itr = r_itr->next;
			}
			
		}

		p_itr = p_itr->next;
		if(p_itr != nullptr){
			r_itr->next = new polynomial();
			r_itr->next->prev = r_itr;			
			r_itr = r_itr->next;
		}
		q_itr = q;
		
	}

	// TODO: Collect like terms
	// optional for now
}

void mult_mono(monomial* a, monomial* b, monomial* result){
	std::vector<bool> found_arr(b->arr.size(), false);
	term* t;

	result->coeff = a->coeff * b->coeff;


	for(int i = 0; i < a->arr.size(); i++){
		t = new term{ .v = a->arr[i]->v, .exp = a->arr[i]->exp };

		for(int j = 0; j < b->arr.size(); j++){
			if(a->arr[i]->v == b->arr[j]->v){
				t->exp += b->arr[j]->exp;
				found_arr[j] = true;
			}
		}

		result->arr.push_back(t);
	}

	for(int i = 0; i < found_arr.size(); i++){
		if(found_arr[i] == false){
			t = new term{ .v = b->arr[i]->v, .exp = b->arr[i]->exp };
			result->arr.push_back(t);
		}
	}
}

double expect_poly(polynomial* p){
	double ret = 0.0;

	polynomial* p_iter = p;
	while (p_iter != nullptr){
		ret += expect_mono(p_iter->m);

		p_iter = p_iter->next;
	}


	return ret;
}

double expect_mono(monomial* m){
	double ret = m->coeff;
	for(int i = 0; i < m->arr.size(); i++){
		if(m->arr[i]->exp%2 != 0){
			ret = 0.0;
			break;
		}
		else{
			ret *= expect_term(m->arr[i]);
		}
	}

	return ret;
}

double expect_term(term* t){
	if(t->v->dist_type == GAUSSIAN){
		if(t->exp % 2 == 0){
			int k = t->exp / 2;
			double numerator = factorial(2*k);

			double denom = std::pow(2.0, k) * factorial(k);
			return numerator / denom;
		}
		else{
			return 0.0;
		}
	}
	else if (t->v->dist_type == UNIFORM) {
		return (std::pow(t->v->b, t->exp + 1) - std::pow(t->v->a, t->exp + 1)) / ((t->v->b - t->v->a) * (t->exp + 1));
	}
}