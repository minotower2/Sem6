#ifndef FUNCTIONS_H
#define FUNCTIONS_H

using func_t = double(*)(double);

func_t get_function(int k);
func_t get_derivative(int k);
const char* get_function_description(int k);

#endif
