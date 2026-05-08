#include "functions.h"
#include "spline_approximation.h"

void make_spline(int n, double* x, double* f, double* c){
	double f_x1_x2 = (f[1] - f[0]) / (x[1] - x[0]);
	c[0] = 0;
	c[1] = 1.;
	c[2] = 0;
	c[3] = f_x1_x2 + (x[0] - x[1]) * ((((f[2] - f[1]) / (x[2] - x[1])) - f_x1_x2) / (x[2] - x[0]));
	//double temp = (f[1] - f[0]) * (x[2] - x[1]) / (x[1] - x[0]);
	for(int i = 1;i < n - 1;i++){
		c[4 * i] = x[i + 1] - x[i];
		c[4 * i + 1] = 2 * (x[i + 1] - x[i - 1]);
		c[4 * i + 2] = x[i] - x[i - 1];
		c[4 * i + 3] = 3 * (((f[i] - f[i - 1]) * (x[i + 1] - x[i]) / (x[i] - x[i - 1])) + (f[i + 1] - f[i]) * (x[i] - x[i - 1]) / (x[i + 1] - x[i]));
		//temp = (f[i + 1] - f[i]) * (x[i] - x[i - 1]) / (x[i + 1] - x[i]);
		//c[4 * i + 3] += temp;
		//c[4 * i + 3] *= 3;
	}
	f_x1_x2 = (f[n - 3] - f[n - 2]) / (x[n - 3] - x[n - 2]);
	c[4 * (n - 1)] = 0;
	c[4 * (n - 1) + 1] = 1.;
	c[4 * (n - 1) + 2] = 0;
	c[4 * (n - 1) + 3] = f_x1_x2 + (2 * x[n - 1] - x[n - 2] - x[n - 3]) * ((((f[n - 1] - f[n - 2]) / (x[n - 1] - x[n - 2])) - f_x1_x2) / (x[n - 1] - x[n - 3]));
	
	double temp;
	for(int i = 0;i < n - 1;i++){// forward elimination
		temp = 1. / c[4 * i + 1];
		//c[4 * i] *= temp;
		//c[4 * i + 1] = 1.;
		c[4 * i + 2] *= temp;
		c[4 * i + 3] *= temp;
		
		temp = c[4 * (i + 1)];
		c[4 * (i + 1) + 3] -= c[4 * i + 3] * temp;
		c[4 * (i + 1) + 1] -= c[4 * i + 2] * temp;
	}
	
	c[4 * (n - 1) + 1] = c[4 * (n - 1) + 3];
	for(int i = n - 2;i >= 0;i--){
		c[4 * i + 1] = c[4 * i + 3] - (c[4 * (i + 1) + 1] * c[4 * i + 2]);
	}
	
	for(int i = 0;i < n - 1;i++){
		c[4 * i] = f[i];
		temp = (f[i + 1] - f[i]) / (x[i + 1] - x[i]);
		c[4 * i + 2] = (3 * temp - 2 * c[4 * i + 1] - c[4 * (i + 1) + 1]) / (x[i + 1] - x[i]);
		c[4 * i + 3] = (c[4 * i + 1] + c[4 * (i + 1) + 1] - 2 * temp) / ((x[i + 1] - x[i]) * (x[i + 1] - x[i]));
	}
}

double calculate_spline_approximation(double x_0, int n, double* x, double* c){
	int l = 0, r = n, s;
	while(l < r){
		s = (l + r) / 2;
		if(x[s] <= x_0){
			l = s + 1;
		} else {
			r = s;
		}
	}
	if(l > 0){
		l--;
	}
	return c[4 * l] + c[4 * l + 1] * (x_0 - x[l]) + c[4 * l + 2] * (x_0 - x[l]) * (x_0 - x[l]) + c[4 * l + 3] * (x_0 - x[l]) * (x_0 - x[l]) * (x_0 - x[l]);
}

double calculate_spline_discrepancy(double x_0, func_t func, int n, double* x, double* c){
	double res = func(x_0) - calculate_spline_approximation(x_0, n, x, c);
	return (res < 0 ? -res : res);
}