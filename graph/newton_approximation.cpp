#include "functions.h"
#include "newton_approximation.h"

void make_Lagrange_polynomial(int n, double* x, double* f, double* d, double* c){
	c[0] = f[0];
	for(int i = 1;i < n;i++){
		c[2 * i - 1] = d[i - 1];
		c[2 * i] = (f[i] - f[i - 1]) / (x[i] - x[i - 1]);
	}
	c[2 * n - 1] = d[n - 1];
	for(int i = 2;i < 2 * n;i++){
		for(int j = 2 * n - 1;j >= i;j--){
			c[j] = (c[j] - c[j - 1]) / (x[j / 2] - x[(j - i) / 2]);
		}
	}
}

double calculate_newton_approximation(double x_0, int n, double* x, double* c){
	double approx_val = c[2 * n - 1];
	for(int i = 2 * n - 2;i >= 0;i--){
		approx_val *= (x_0 - x[i / 2]);
		approx_val += c[i];
	}
	return approx_val;
}

double calculate_newton_discrepancy(double x_0, func_t func, int n, double* x, double* c){
	double res = func(x_0) - calculate_newton_approximation(x_0, n, x, c);
	return (res < 0 ? -res : res);
}