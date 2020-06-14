
#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include "lp_lib.h"
using namespace std;

const double Error = 1e-10;
int number_of_cash_flows;
vector <double> price_list;
vector <int> maturity_list;
vector <double> yield_to_maturity;
vector <double> duration;
vector <double> convexity;
double debt_obligation_amount;
double time_when_debt_is_due;
vector <double> percentage_of_cash_flow_to_meet_debt_obligation;
vector <vector <double> > cash_flow;

double function(vector <double> cash_flow, double price, int maturity, double rate)
{
	//computes f(r).
	double f = price * pow(rate + 1, maturity);
	for (int i = 1; i <= maturity; i++)
	{
		f = f - cash_flow[i-1] * pow(rate + 1, maturity - i);
	}
	return f;
}

double derivative_function(vector <double> cash_flow, double price, int maturity, double rate)
{
	//computes f'(r).
	double f_1 = maturity * price * pow(rate + 1, maturity - 1);
	for (int i = 1; i <= maturity - 1; i++)
	{
		f_1 = f_1 - cash_flow[i - 1] * (maturity - i) * pow(rate + 1, maturity - i - 1);
	}
	return f_1;
}

double Newton_Raphson(vector <double> cash_flow, double price, int maturity, double rate)
{
	//using Newton - Raphson method to find the only root(the YTM of each bond or cash flow).
	double r = rate;
	double e = 1; //e is the difference between input r each time and the root.
	while(e > Error) //when e < Error, we think the r is approximately the same as root.
	{
		r = r - function(cash_flow, price, maturity, r) / derivative_function(cash_flow, price, maturity, r);
		e = function(cash_flow, price, maturity, r) / derivative_function(cash_flow, price, maturity, r);
	}
	return r;
}

double get_duration(vector <double> cash_flow, double price, int maturity, double rate)
{
	//computes duration of each bond(cash flow).
	double N = 0;
	for (int i = 1; i <= maturity; i++)
	{
		N = N + ((i * cash_flow[i - 1]) / pow(1 + rate, i));
	}
	N = N / price;
	return N;
}

double get_convexity(vector <double> cash_flow, double price, int maturity, double rate)
{
	//computes convexity of each bond(cash flow).
	double convexity = 0;
	for (int i = 1; i <= maturity; i++)
	{
		convexity = convexity + ((i * (i+1) * cash_flow[i - 1]) / pow(1 + rate, i + 2));
	}
	convexity = convexity / price;
	return convexity;
}

double present_value_of_debt()
{
	//computes the present value of debt obligation using the average value of the YTMs.
	double sum = 0;
	double aver_rate = 0;
	double pre_debt = 0;
	for (int i = 0; i < number_of_cash_flows; i++)
	{
		sum = sum + yield_to_maturity[i];
	}
	aver_rate = sum / number_of_cash_flows;
	pre_debt = debt_obligation_amount / pow(1 + aver_rate, time_when_debt_is_due);
	return pre_debt;
}

void print_data(char *filename)
{
	cout << "Input File: " << filename << endl;
	cout << "We owe " << debt_obligation_amount << " in " << time_when_debt_is_due << " years" << endl;
	cout << "Number of Cash Flows: " << number_of_cash_flows << endl;
	for (int i = 0; i < number_of_cash_flows; i++) 
	{
		cout << "-------------------------------------------------------------------" << endl;
		cout << "Cash Flow #" << i+1 << endl;
		cout << "Price = " << price_list[i] << endl;
		cout << "Maturity = " << maturity_list[i] << endl;
		cout << "Yield to Maturity = " << yield_to_maturity[i] << endl;
		cout << "Duration = " << duration[i] << endl;
		cout << "Convexity = " << convexity[i] << endl;
		cout << "Percentage of Face Value that would meet the obligation = " << percentage_of_cash_flow_to_meet_debt_obligation[i] << endl;
	}
	cout << "-------------------------------------------------------------------" << endl;
}

void get_data(char* argv[])
{
	double element; //set a double type variable to transport value.
	int ele; //set a int type variable to transport value.

	//reads the data from the file identified as follow.
	ifstream input_filename(argv[1]);
	input_filename >> number_of_cash_flows;
	int q = 0;
	for (int j = 1; j <= number_of_cash_flows; j++)
	{
		input_filename >> element;
		price_list.push_back(element);
		input_filename >> ele;
		maturity_list.push_back(ele);
		cash_flow.push_back(vector <double>());
		for (int i = 1; i <= ele; i++)
		{
			input_filename >> element;
			cash_flow[j - 1].push_back(element);
		}
	}
	input_filename >> debt_obligation_amount;
	input_filename >> time_when_debt_is_due;

	//compute the bonds'(cash flows') YTMs, Durations, Convexities and Percentage of their Face Values that would meet the obligation as follow.
	double yield = 0;
	double pre = 0;
	for (int k = 1; k <= number_of_cash_flows; k++)
	{
		yield = Newton_Raphson(cash_flow[k - 1], price_list[k - 1], maturity_list[k - 1], 0.06);
		yield_to_maturity.push_back(yield);
		duration.push_back(get_duration(cash_flow[k - 1], price_list[k - 1], maturity_list[k - 1], yield));
		convexity.push_back(get_convexity(cash_flow[k - 1], price_list[k - 1], maturity_list[k - 1], yield));
	}
	pre = present_value_of_debt();
	for (int l = 0; l < number_of_cash_flows; l++)
	{
		percentage_of_cash_flow_to_meet_debt_obligation.push_back(pre / price_list[l]);
	}
}

void get_optimal_portfolio()
{
	lprec *lp;          //set lp model.
	double* cost;       //since the default of lpsolve is to minimize the target formula, so in order to maximize it, we set a pointer to save opposite number of convexities.
	cost = new double [number_of_cash_flows];
	double* solution;   //set a pointer to put the solutions of lp.
	solution = new double[number_of_cash_flows];
	for (int i = 0; i < number_of_cash_flows; i++)   //assign the values to cost.
	{
		cost[i] = (-1) * convexity[i];
	}
	lp = make_lp(0, number_of_cash_flows);   //create a model with 0 rows and the number of cash flows columns(variables).
	set_verbose(lp, 3);                      //only show something important on the screen.
	{
		//create the first row as a constraint: Make the sum of solutions equals to 1, since necessary bonds(cash flows) are parts of the portfolio.
		double* row;
		row = new double[number_of_cash_flows];
		for (int i = 1; i <= number_of_cash_flows; i++)
		{
			row[i] = 1;
		}
		add_constraint(lp, row, EQ, 1);
	}
	{
		//create the second row as a constraint: Make the portfolio's duration equals to the debt obligation's.
		double* row;
		row = new double[number_of_cash_flows];
		for (int i = 1; i <= number_of_cash_flows; i++)
		{
			row[i] = duration[i-1];
		}
		add_constraint(lp, row, EQ, time_when_debt_is_due);
	}
	{
		//set the problem row. Our target is to maximize the portfolio's convexity.
		double* row;
		row = new double[number_of_cash_flows];
		for (int i = 1; i <= number_of_cash_flows; i++)
		{
			row[i] = cost[i-1];
		}
		set_obj_fn(lp, row);
	}
	print_lp(lp); //output the basic information about our lp model.
	int ret;
	ret = solve(lp);
	if (ret == 0) //if the model has solutions, then output them.
	{
		cout << "Largest Convexity we can get is: " << (-1) * get_objective(lp) << endl;   //return the value of the objective of the last solve(the largest convexity in possible portfolios).
		cout << "Optimal portfolio:" << endl;
		get_variables(lp, solution);                   //put the solutions into a well-setting pointer to contain the values.
		for (int i = 0; i < number_of_cash_flows; i++) //show the solutions.
		{
			cout << "%Cash Flow:" << i + 1 << "  " << solution[i] << endl;
		}
		cout << "-------------------------------------------------------------------" << endl;
		cout << "To immunize against small changes in 'r' for each $1 of PV, you should buy" << endl;
		for (int i = 0; i < number_of_cash_flows; i++) //show the final composition of portfolio with the best immunization.
		{
			if (solution[i] != 0)
				cout << "$" << solution[i] << " of Cash Flow#" << i + 1 << endl;
		}
		cout << "If you need to immunize for a larger PV-value, just buy an appropriate proportion" << endl;
		cout << "-------------------------------------------------------------------" << endl;
		//show the value of each bond(cash flow) that you should buy for various PV's immunization.
		cout << "For example, if you want to immunize for $500 of PV, buy" << endl;
		for (int i = 0; i < number_of_cash_flows; i++)
		{
			if (solution[i] != 0)
				cout << "$" << solution[i] * 500 << " of Cash Flow#" << i + 1 << endl;
		}
		cout << "-------------------------------------------------------------------" << endl;
		cout << "For example, if you want to immunize for $750 of PV, buy" << endl;
		for (int i = 0; i < number_of_cash_flows; i++)
		{
			if (solution[i] != 0)
				cout << "$" << solution[i] * 750 << " of Cash Flow#" << i + 1 << endl;
		}
		cout << "-------------------------------------------------------------------" << endl;
		cout << "For example, if you want to immunize for $1000 of PV, buy" << endl;
		for (int i = 0; i < number_of_cash_flows; i++)
		{
			if (solution[i] != 0)
				cout << "$" << solution[i] * 1000 << " of Cash Flow#" << i + 1 << endl;
		}
		cout << "-------------------------------------------------------------------" << endl;
		cout << "For example, if you want to immunize for $1009.36 of PV, buy" << endl;
		for (int i = 0; i < number_of_cash_flows; i++)
		{
			if (solution[i] != 0)
				cout << "$" << solution[i] * 1009.36 << " of Cash Flow#" << i + 1 << endl;
		}
		cout << "-------------------------------------------------------------------" << endl;
	}
	else
	{
		//when there is no solution of lp model.
		cout << "There is no potfolio that meets the duration constraint of " << time_when_debt_is_due << " years." << endl;
	}
}
	
int main (int argc, char* argv[])
{
	if (argc == 1) {
		cout << "Input filename missing" << endl;
	}
	else 
	{
		get_data(argv);
		print_data(argv[1]);
		get_optimal_portfolio();
	}
	return (0);
}

