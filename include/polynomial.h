#ifndef _POLY_UTILS_H
#define _POLY_UTILS_H

#include <math.h>
#include <vector>
#include <iostream>

#include "../include/enums.h"

class var {
	public:
	int id;
	double b;
	double a;

	public:
	var(double b, double a, int id);
};

class term {
	public:
	var* v;
	int exp;
};

class monomial{
	public:
	double coeff;
	std::vector<term*> arr;

	void print();
	monomial();
	~monomial();
};

class polynomial;

class polynomial{
	public:
	polynomial();
	monomial* m;
	polynomial* next;
	polynomial* prev;
	int max_exp;

	void print();
	polynomial* copy();
	bool equals(polynomial* p);
	~polynomial();
};

void faster_tensor_prod(std::vector<std::vector<polynomial*>>& univariate_set, std::vector<polynomial*>& result, int max_order);
void tensor_prod(std::vector<std::vector<polynomial*>>& univariate_set, std::vector<polynomial*>& result, int max_order);
void tensor_prod(std::vector<polynomial*>& set_a, std::vector<polynomial*>& set_b, std::vector<polynomial*>& result, int max_order);

void mult_poly(polynomial* p, polynomial* q, polynomial* result);
void mult_mono(monomial* a, monomial* b, monomial* result);

double expect_poly(polynomial* p, RandVarDist dist_type);
double expect_mono(monomial* m, RandVarDist dist_type);
double expect_term(term* t, RandVarDist dist_type);

#endif
