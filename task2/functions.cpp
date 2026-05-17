#include <cmath>
#include "functions.h"

double f_0(double x){
	(void)x;
	return 1;
}

double f_1(double x){
	return x;
}

double f_2(double x){
	return x * x;
}

double f_3(double x){
	return x * x * x;
}

double f_4(double x){
	return x * x * x * x;
}

double f_5(double x){
	return std::exp(x);
}

double f_6(double x){
	return 1. / (25 * x * x + 1);
}

/* ------------- */

double d_0(double x){
	(void)x;
	return 0;
}

double d_1(double x){
	(void)x;
	return 1;
}

double d_2(double x){
	return 2 * x;
}

double d_3(double x){
	return 3 * x * x;
}

double d_4(double x){
	return 4 * x * x * x;
}

double d_5(double x){
	return std::exp(x);
}

double d_6(double x){
	return 0 - (50 * x) / (625 * x * x * x * x + 50 * x * x + 1);
}

/* ------------- */

func_t get_function(int k){
	func_t f[] = {f_0, f_1, f_2, f_3, f_4, f_5, f_6};
	if(k < 0 || static_cast<unsigned int>(k) >= (sizeof(f) / sizeof(f[0]))){
		k = 0;
	}
	return f[k];
}

func_t get_derivative(int k){
	func_t f[] = {d_0, d_1, d_2, d_3, d_4, d_5, d_6};
	if(k < 0 || static_cast<unsigned int>(k) >= (sizeof(f) / sizeof(f[0]))){
		k = 0;
	}
	return f[k];
}

const char* get_function_description(int k){
	switch(k){
		default: return "f(x) = 1";// k == 0
		case 1: return "f(x) = x";
		case 2: return "f(x) = x^2";
		case 3: return "f(x) = x^3";
		case 4: return "f(x) = x^4";
		case 5: return "f(x) = e^x";
		case 6: return "f(x) = 1/(25x^2 + 1)";
	}
	return "";
}