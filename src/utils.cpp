#include "../include/utils.h"

double factorial(int n){
	double retval = 1.0;

	if(n != 0){
		for(int i = 1; i <= n; i++){
			retval *= i;
		}
	}
	
	return retval;
}