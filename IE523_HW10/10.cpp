#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include "lp_lib.h"
using namespace std;

double N(const double& z) {
	if (z > 6.0) { return 1.0; }; // this guards against overflow 
	if (z < -6.0) { return 0.0; };
	double b1 = 0.31938153;
	double b2 = -0.356563782;
	double b3 = 1.781477937;
	double b4 = -1.821255978;
	double b5 = 1.330274429;
	double p = 0.2316419;
	double c2 = 0.3989423;
	double a = fabs(z);
	double t = 1.0 / (1.0 + a * p);
	double b = c2 * exp((-z) * (z / 2.0));
	double n = ((((b5 * t + b4) * t + b3) * t + b2) * t + b1) * t;
	n = 1.0 - b * n;
	if (z < 0.0) n = 1.0 - n;
	return n;
};

double up_factor, risk_free_rate, strike_price, R;
double initial_stock_price, expiration_time, volatility;
int** call_memoization;
int** put_memoization;
int no_of_divisions;
lprec* lp_call, * lp_put, * lp_call_via_memoization, * lp_put_via_memoization;

double MAX(double a, double b)
{
	return ((b < a) ? a : b);
}

double option_price_put_black_scholes(const double& S,      // spot price
	const double& K,      
	const double& r,      
	const double& sigma,  // volatility
	const double& time)   // time to maturity 
{
	double time_sqrt = sqrt(time);
	double d1 = (log(S / K) + r * time) / (sigma * time_sqrt) + 0.5 * sigma * time_sqrt;
	double d2 = d1 - (sigma * time_sqrt);
	return K * exp(-r * time) * N(-d2) - S * N(-d1);
};

double option_price_call_black_scholes(const double& S,       // spot (underlying) price
	const double& K,      
	const double& r,       
	const double& sigma,   
	const double& time) 
{  
	double time_sqrt = sqrt(time);
	double d1 = (log(S / K) + r * time) / (sigma * time_sqrt) + 0.5 * sigma * time_sqrt;
	double d2 = d1 - (sigma * time_sqrt);
	return S * N(d1) - K * exp(-r * time) * N(d2);
};

void create_LP_for_european_call_option(int k, int i)
{
	if (k == no_of_divisions - 1)
	{
		double* row1 = new double[no_of_divisions * (no_of_divisions + 1) + 1];
		double* row2 = new double[no_of_divisions * (no_of_divisions + 1) + 1];
		for (int j = 0; j < no_of_divisions * (no_of_divisions + 1) + 1; j++) 
		{
			row1[j] = 0.0;
			row2[j] = 0.0;
		}

		row1[(k * k) + k + (i + k) + 1] = initial_stock_price * pow(up_factor, ((double)i + 1));
		row1[(k * k) + k + (i + k) + 2] = -R;
		row2[(k * k) + k + (i + k) + 1] = initial_stock_price * pow(up_factor, ((double)i - 1));
		row2[(k * k) + k + (i + k) + 2] = -R;

		add_constraint(lp_call, row1, GE, MAX(0.0, (initial_stock_price * pow(up_factor, ((double)i + 1))) - strike_price));
		add_constraint(lp_call, row2, GE, MAX(0.0, (initial_stock_price * pow(up_factor, ((double)i - 1))) - strike_price));
	}
	else 
	{
		double* row1 = new double[no_of_divisions * (no_of_divisions + 1) + 1];
		double* row2 = new double[no_of_divisions * (no_of_divisions + 1) + 1];
		for (int j = 0; j < no_of_divisions * (no_of_divisions + 1) + 1; j++) 
		{
			row1[j] = 0.0;
			row2[j] = 0.0;
		}

		row1[(k * k) + k + (i + k) + 1] = initial_stock_price * pow(up_factor, ((double)i + 1));
		row1[(k * k) + k + (i + k) + 2] = -R;
		row2[(k * k) + k + (i + k) + 1] = initial_stock_price * pow(up_factor, ((double)i - 1));
		row2[(k * k) + k + (i + k) + 2] = -R;

		row1[((k + 1) * (k + 1)) + (k + 1) + (i + 1 + k + 1) + 1] = -initial_stock_price * pow(up_factor, ((double)i + 1));
		row1[((k + 1) * (k + 1)) + (k + 1) + (i + 1 + k + 1) + 2] = 1;
		row2[((k + 1) * (k + 1)) + (k + 1) + (i - 1 + k + 1) + 1] = -initial_stock_price * pow(up_factor, ((double)i - 1));
		row2[((k + 1) * (k + 1)) + (k + 1) + (i - 1 + k + 1) + 2] = 1;

		add_constraint(lp_call, row1, GE, 0);
		add_constraint(lp_call, row2, GE, 0);

		create_LP_for_european_call_option(k + 1, i + 1);
		create_LP_for_european_call_option(k + 1, i - 1);
	}
}

void create_LP_for_european_call_option_via_memoization(int k, int i)
{
	if (-1 == call_memoization[k][(i + k) / 2])
	{
		if (k == no_of_divisions - 1)
		{
			double* row1 = new double[no_of_divisions * (no_of_divisions + 1) + 1];
			double* row2 = new double[no_of_divisions * (no_of_divisions + 1) + 1];
			for (int j = 0; j < no_of_divisions * (no_of_divisions + 1) + 1; j++) 
			{
				row1[j] = 0.0;
				row2[j] = 0.0;
			}

			row1[(k * k) + k + (i + k) + 1] = initial_stock_price * pow(up_factor, ((double)i + 1));
			row1[(k * k) + k + (i + k) + 2] = -R;
			row2[(k * k) + k + (i + k) + 1] = initial_stock_price * pow(up_factor, ((double)i - 1));
			row2[(k * k) + k + (i + k) + 2] = -R;

			add_constraint(lp_call_via_memoization, row1, GE, MAX(0.0, (initial_stock_price * pow(up_factor, ((double)i + 1))) - strike_price));
			add_constraint(lp_call_via_memoization, row2, GE, MAX(0.0, (initial_stock_price * pow(up_factor, ((double)i - 1))) - strike_price));
		}
		else 
		{
			double* row1 = new double[no_of_divisions * (no_of_divisions + 1) + 1];
			double* row2 = new double[no_of_divisions * (no_of_divisions + 1) + 1];
			for (int j = 0; j < no_of_divisions * (no_of_divisions + 1) + 1; j++) 
			{
				row1[j] = 0.0;
				row2[j] = 0.0;
			}

			row1[(k * k) + k + (i + k) + 1] = initial_stock_price * pow(up_factor, ((double)i + 1));
			row1[(k * k) + k + (i + k) + 2] = -R;
			row2[(k * k) + k + (i + k) + 1] = initial_stock_price * pow(up_factor, ((double)i - 1));
			row2[(k * k) + k + (i + k) + 2] = -R;

			row1[((k + 1) * (k + 1)) + (k + 1) + (i + 1 + k + 1) + 1] = -initial_stock_price * pow(up_factor, ((double)i + 1));
			row1[((k + 1) * (k + 1)) + (k + 1) + (i + 1 + k + 1) + 2] = 1;
			row2[((k + 1) * (k + 1)) + (k + 1) + (i - 1 + k + 1) + 1] = -initial_stock_price * pow(up_factor, ((double)i - 1));
			row2[((k + 1) * (k + 1)) + (k + 1) + (i - 1 + k + 1) + 2] = 1;

			add_constraint(lp_call_via_memoization, row1, GE, 0);
			add_constraint(lp_call_via_memoization, row2, GE, 0);

			create_LP_for_european_call_option_via_memoization(k + 1, i + 1);
			create_LP_for_european_call_option_via_memoization(k + 1, i - 1);
		}

		call_memoization[k][(i + k) / 2] = 1;
	}
}

void create_LP_for_european_put_option(int k, int i)
{
	if (k == no_of_divisions - 1)
	{
		// creating the constraint row & initializing with all zeros
		double* row1 = new double[no_of_divisions * (no_of_divisions + 1) + 1];
		double* row2 = new double[no_of_divisions * (no_of_divisions + 1) + 1];
		for (int j = 0; j < no_of_divisions * (no_of_divisions + 1) + 1; j++) 
		{
			row1[j] = 0.0;
			row2[j] = 0.0;
		}

		row1[(k * k) + k + (i + k) + 1] = -initial_stock_price * pow(up_factor, ((double)i + 1));
		row1[(k * k) + k + (i + k) + 2] = R;
		row2[(k * k) + k + (i + k) + 1] = -initial_stock_price * pow(up_factor, ((double)i - 1));
		row2[(k * k) + k + (i + k) + 2] = R;

		add_constraint(lp_put, row1, GE, MAX(0.0, (strike_price - initial_stock_price * pow(up_factor, ((double)i + 1)))));
		add_constraint(lp_put, row2, GE, MAX(0.0, (strike_price - initial_stock_price * pow(up_factor, ((double)i - 1)))));
	}
	else
	{
		double* row1 = new double[no_of_divisions * (no_of_divisions + 1) + 1];
		double* row2 = new double[no_of_divisions * (no_of_divisions + 1) + 1];
		for (int j = 0; j < no_of_divisions * (no_of_divisions + 1) + 1; j++)
		{
			row1[j] = 0.0;
			row2[j] = 0.0;
		}

		row1[(k * k) + k + (i + k) + 1] = -initial_stock_price * pow(up_factor, ((double)i + 1));
		row1[(k * k) + k + (i + k) + 2] = R;
		row2[(k * k) + k + (i + k) + 1] = -initial_stock_price * pow(up_factor, ((double)i - 1));
		row2[(k * k) + k + (i + k) + 2] = R;

		row1[((k + 1) * (k + 1)) + (k + 1) + (i + 1 + k + 1) + 1] = initial_stock_price * pow(up_factor, ((double)i + 1));
		row1[((k + 1) * (k + 1)) + (k + 1) + (i + 1 + k + 1) + 2] = -1;
		row2[((k + 1) * (k + 1)) + (k + 1) + (i - 1 + k + 1) + 1] = initial_stock_price * pow(up_factor, ((double)i - 1));
		row2[((k + 1) * (k + 1)) + (k + 1) + (i - 1 + k + 1) + 2] = -1;

		add_constraint(lp_put, row1, GE, 0);
		add_constraint(lp_put, row2, GE, 0);

		create_LP_for_european_put_option(k + 1, i + 1);
		create_LP_for_european_put_option(k + 1, i - 1);
	}
}

void create_LP_for_european_put_option_via_memoization(int k, int i)
{
	if (-1 == put_memoization[k][(i + k) / 2])
	{
		if (k == no_of_divisions - 1)
		{
			double* row1 = new double[no_of_divisions * (no_of_divisions + 1) + 1];
			double* row2 = new double[no_of_divisions * (no_of_divisions + 1) + 1];
			for (int j = 0; j < no_of_divisions * (no_of_divisions + 1) + 1; j++) 
			{
				row1[j] = 0.0;
				row2[j] = 0.0;
			}

			row1[(k * k) + k + (i + k) + 1] = -initial_stock_price * pow(up_factor, ((double)i + 1));
			row1[(k * k) + k + (i + k) + 2] = R;
			row2[(k * k) + k + (i + k) + 1] = -initial_stock_price * pow(up_factor, ((double)i - 1));
			row2[(k * k) + k + (i + k) + 2] = R;

			add_constraint(lp_put_via_memoization, row1, GE, MAX(0.0, (strike_price - initial_stock_price * pow(up_factor, ((double)i + 1)))));
			add_constraint(lp_put_via_memoization, row2, GE, MAX(0.0, (strike_price - initial_stock_price * pow(up_factor, ((double)i - 1)))));
		}
		else
		{
			double* row1 = new double[no_of_divisions * (no_of_divisions + 1) + 1];
			double* row2 = new double[no_of_divisions * (no_of_divisions + 1) + 1];
			for (int j = 0; j < no_of_divisions * (no_of_divisions + 1) + 1; j++)
			{
				row1[j] = 0.0;
				row2[j] = 0.0;
			}

			row1[(k * k) + k + (i + k) + 1] = -initial_stock_price * pow(up_factor, ((double)i + 1));
			row1[(k * k) + k + (i + k) + 2] = R;
			row2[(k * k) + k + (i + k) + 1] = -initial_stock_price * pow(up_factor, ((double)i - 1));
			row2[(k * k) + k + (i + k) + 2] = R;

			row1[((k + 1) * (k + 1)) + (k + 1) + (i + 1 + k + 1) + 1] = initial_stock_price * pow(up_factor, ((double)i + 1));
			row1[((k + 1) * (k + 1)) + (k + 1) + (i + 1 + k + 1) + 2] = -1;
			row2[((k + 1) * (k + 1)) + (k + 1) + (i - 1 + k + 1) + 1] = initial_stock_price * pow(up_factor, ((double)i - 1));
			row2[((k + 1) * (k + 1)) + (k + 1) + (i - 1 + k + 1) + 2] = -1;

			add_constraint(lp_put_via_memoization, row1, GE, 0);
			add_constraint(lp_put_via_memoization, row2, GE, 0);

			create_LP_for_european_put_option_via_memoization(k + 1, i + 1);
			create_LP_for_european_put_option_via_memoization(k + 1, i - 1);
		}

		put_memoization[k][(i + k) / 2] = 1;
	}
}

void set_up_and_solve_the_LP_for_the_call_option()
{
	lp_call = make_lp(0, no_of_divisions * (no_of_divisions + 1));
	set_verbose(lp_call, 3);
	create_LP_for_european_call_option(0, 0);
	{
		double* row = new double[no_of_divisions * (no_of_divisions + 1) + 1];
		for (int i = 0; i < no_of_divisions * (no_of_divisions + 1) + 1; i++)
			row[i] = 0.0;
		row[1] = initial_stock_price;
		row[2] = -1;
		set_obj_fn(lp_call, row);
	}
	{
		int error_code = solve(lp_call);
		if (0 == error_code)
			cout << "Call Price according the LP formulation = " << get_objective(lp_call) << endl;
		else 
		{
			cout << "LP solve ran into problems; Error Code = " << error_code << endl;
			cout << "Look up http://lpsolve.sourceforge.net/5.1/solve.htm for error code details" << endl;
		}
	}
	delete_lp(lp_call);
}

void set_up_and_solve_the_LP_for_the_call_option_via_memoization()
{
	lp_call_via_memoization = make_lp(0, no_of_divisions * (no_of_divisions + 1));
	set_verbose(lp_call_via_memoization, 3);
	create_LP_for_european_call_option_via_memoization(0, 0);
	{
		double* row = new double[no_of_divisions * (no_of_divisions + 1) + 1];
		for (int i = 0; i < no_of_divisions * (no_of_divisions + 1) + 1; i++)
			row[i] = 0.0;
		row[1] = initial_stock_price;
		row[2] = -1;
		set_obj_fn(lp_call_via_memoization, row);
	}
	{
		int error_code = solve(lp_call_via_memoization);
		if (0 == error_code)
			cout << "Call Price according the LP formulation = " << get_objective(lp_call_via_memoization) << endl;
		else 
		{
			cout << "LP solve ran into problems; Error Code = " << error_code << endl;
			cout << "Look up http://lpsolve.sourceforge.net/5.1/solve.htm for error code details" << endl;
		}
	}
	delete_lp(lp_call_via_memoization);
}

void set_up_and_solve_the_LP_for_the_put_option()
{
	lp_put = make_lp(0, no_of_divisions * (no_of_divisions + 1));
	set_verbose(lp_put, 3);
	create_LP_for_european_put_option(0, 0);
	{
		double* row = new double[no_of_divisions * (no_of_divisions + 1) + 1];
		for (int i = 0; i < no_of_divisions * (no_of_divisions + 1) + 1; i++)
			row[i] = 0.0;
		row[1] = -initial_stock_price;
		row[2] = 1;
		set_obj_fn(lp_put, row);
	}
	{
		int error_code = solve(lp_put);
		if (0 == error_code)
			cout << "Put Price according the LP formulation = " << get_objective(lp_put) << endl;
		else 
		{
			cout << "LP solve ran into problems; Error Code = " << error_code << endl;
			cout << "Look up http://lpsolve.sourceforge.net/5.1/solve.htm for error code details" << endl;
		}
	}
	delete_lp(lp_put);
}

void set_up_and_solve_the_LP_for_the_put_option_via_memoization()
{
	lp_put_via_memoization = make_lp(0, no_of_divisions * (no_of_divisions + 1));
	set_verbose(lp_put_via_memoization, 3);
	create_LP_for_european_put_option_via_memoization(0, 0);
	{
		double* row = new double[no_of_divisions * (no_of_divisions + 1) + 1];
		for (int i = 0; i < no_of_divisions * (no_of_divisions + 1) + 1; i++)
			row[i] = 0.0;
		row[1] = -initial_stock_price;
		row[2] = 1;
		set_obj_fn(lp_put_via_memoization, row);
	}
	{
		int error_code = solve(lp_put_via_memoization);
		if (0 == error_code)
			cout << "Put Price according the LP formulation = " << get_objective(lp_put_via_memoization) << endl;
		else 
		{
			cout << "LP solve ran into problems; Error Code = " << error_code << endl;
			cout << "Look up http://lpsolve.sourceforge.net/5.1/solve.htm for error code details" << endl;
		}
	}
	delete_lp(lp_put_via_memoization);
}


int main(int argc, char* argv[])
{
	sscanf_s(argv[1], "%lf", &expiration_time);
	sscanf_s(argv[2], "%d", &no_of_divisions);
	sscanf_s(argv[3], "%lf", &risk_free_rate);
	sscanf_s(argv[4], "%lf", &volatility);
	sscanf_s(argv[5], "%lf", &initial_stock_price);
	sscanf_s(argv[6], "%lf", &strike_price);

	call_memoization = new int* [no_of_divisions];
	for (int i = 0; i < no_of_divisions; i++)
		call_memoization[i] = new int[no_of_divisions];
	put_memoization = new int* [no_of_divisions];
	for (int i = 0; i < no_of_divisions; i++)
		put_memoization[i] = new int[no_of_divisions];
	for (int i = 0; i < no_of_divisions; i++)
	{
		for (int j = 0; j < no_of_divisions; j++)
		{
			call_memoization[i][j] = -1;
			put_memoization[i][j] = -1;
		}
	}

	up_factor = exp(volatility * sqrt(expiration_time / ((double)no_of_divisions)));
	R = exp(risk_free_rate * expiration_time / ((double)no_of_divisions));

	cout << "--------------------------------------" << endl;
	cout << "European Option Pricing via Linear Programming with and without memoization" << endl;
	cout << "Expiration Time (Years) = " << expiration_time << endl;
	cout << "Number of Divisions = " << no_of_divisions << endl;
	cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
	cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
	cout << "Initial Stock Price = " << initial_stock_price << endl;
	cout << "Strike Price = " << strike_price << endl;
	cout << "R = " << R << endl;
	cout << "Up-Factor = " << up_factor << endl;
	cout << "--------------------------------------" << endl;
	cout << "Pricing the Call with Memoization" << endl;
	clock_t call_memo_start = clock();
	set_up_and_solve_the_LP_for_the_call_option_via_memoization();
	clock_t call_memo_end = clock();
	cout << "It took " << (long double)(call_memo_end - call_memo_start) / CLOCKS_PER_SEC << " seconds to complete" << endl;
	cout << "--------------------------------------" << endl;
	cout << "Pricing the Call without Memoization" << endl;
	clock_t call_start = clock();
	set_up_and_solve_the_LP_for_the_call_option();
	clock_t call_end = clock();
	cout << "It took " << (long double)(call_end - call_start) / CLOCKS_PER_SEC << " seconds to complete" << endl;
	cout << "--------------------------------------" << endl;
	cout << "Call Price according to Black-Scholes = " <<
		option_price_call_black_scholes(initial_stock_price, strike_price, risk_free_rate, volatility, expiration_time) << endl;
	cout << "--------------------------------------" << endl;
	cout << "                                      " << endl;
	cout << "--------------------------------------" << endl;
	cout << "Pricing the Put with Memoization" << endl;
	clock_t put_memo_start = clock();
	set_up_and_solve_the_LP_for_the_put_option_via_memoization();
	clock_t put_memo_end = clock();
	cout << "It took " << (long double)(put_memo_end - put_memo_start) / CLOCKS_PER_SEC << " seconds to complete" << endl;
	cout << "--------------------------------------" << endl;
	cout << "Pricing the Put without Memoization" << endl;
	clock_t put_start = clock();
	set_up_and_solve_the_LP_for_the_put_option();
	clock_t put_end = clock();
	cout << "It took " << (long double)(put_end - put_start) / CLOCKS_PER_SEC << " seconds to complete" << endl;
	cout << "--------------------------------------" << endl;
	cout << "Put Price according to Black-Scholes = " <<
		option_price_put_black_scholes(initial_stock_price, strike_price, risk_free_rate, volatility, expiration_time) << endl;
	cout << "--------------------------------------" << endl;
}