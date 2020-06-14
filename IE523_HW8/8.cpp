
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <time.h>
#include "newmat.h"
using namespace std;

double up_factor, uptick_prob, downtick_prob, risk_free_rate, strike_price;
double initial_stock_price, expiration_time, volatility, R;
int no_of_divisions;

//Part 1: Pricing an American-Option using a Trinomial Model by Recursion
double american_call_option(int k, int i, double current_stock_price)
{
	if (k == no_of_divisions)
		return max(0.0, (current_stock_price - strike_price));
	else
		return max((current_stock_price - strike_price),
		(uptick_prob * american_call_option(k + 1, i + 1, current_stock_price * up_factor) +
			(1- uptick_prob-downtick_prob)* american_call_option(k + 1, i, current_stock_price)+
			downtick_prob * american_call_option(k + 1, i - 1, current_stock_price / up_factor)) / R);
}

double american_put_option(int k, int i, double current_stock_price)
{
	if (k == no_of_divisions)
		return max(0.0, (strike_price - current_stock_price));
	else
		return max((strike_price - current_stock_price),
		(uptick_prob * american_put_option(k + 1, i + 1, current_stock_price * up_factor) +
			(1 - uptick_prob - downtick_prob) * american_put_option(k + 1, i, current_stock_price) +
			downtick_prob * american_put_option(k + 1, i - 1, current_stock_price / up_factor)) / R);
}

//Part 2: Pricing an American-Option using a Trinomial Model by Dynamic Programming
double max(double a, double b) {
	return (b < a) ? a : b;
}

double american_call_option_dyn_prog()
{
	// create the probability matrix
	Matrix transition_probability(2 * no_of_divisions + 1, 2 * no_of_divisions + 1);

	for (int i = 1; i <= 2 * no_of_divisions + 1; i++)
		for (int j = 1; j <= 2 * no_of_divisions + 1; j++)
			transition_probability(i, j) = 0.0;

	// boundary values of the probabilities need to be entered
	transition_probability(1, 1) = 1.0 - uptick_prob;
	transition_probability(1, 2) = uptick_prob;
	transition_probability(2 * no_of_divisions + 1, 2 * no_of_divisions + 1) = 1.0 - downtick_prob;
	transition_probability(2 * no_of_divisions + 1, 2 * no_of_divisions) = downtick_prob;

	for (int i = 2; i <= 2 * no_of_divisions; i++)
		for (int j = 1; j <= 2 * no_of_divisions + 1; j++) 
		{
			if (j == (i - 1))
				transition_probability(i, j) = downtick_prob;
			if (j == i)
				transition_probability(i, j) = 1.0 - uptick_prob - downtick_prob;
			if (j == (i + 1))
				transition_probability(i, j) = uptick_prob;
		}

	Matrix V_t(2 * no_of_divisions + 1, 1);
	// value at expiration
	for (int i = 1; i <= 2 * no_of_divisions + 1; i++)
		V_t(i, 1) = max(0.0, (initial_stock_price * pow(up_factor, i - 1 - no_of_divisions)) - strike_price);

	// value at intermediate stages
	Matrix discounted_one_step_forward_value(2 * no_of_divisions + 1, 1);
	Matrix value_if_option_is_exercised_now(2 * no_of_divisions + 1, 1);
	for (int i = no_of_divisions; i > 0; i--) 
	{
		// going backwards from expiration to zero-time
		discounted_one_step_forward_value = (transition_probability * V_t) / R;
		for (int j = 1; j <= 2 * no_of_divisions + 1; j++)
			value_if_option_is_exercised_now(j, 1) = max(0, (initial_stock_price * pow(up_factor, j - 1 - no_of_divisions) - strike_price));
		for (int j = 1; j <= 2 * no_of_divisions + 1; j++)
			V_t(j, 1) = max(value_if_option_is_exercised_now(j, 1), discounted_one_step_forward_value(j, 1));
	}
	return (V_t(no_of_divisions + 1, 1));
}

double american_put_option_dyn_prog()
{
	// create the probability matrix
	Matrix transition_probability(2 * no_of_divisions + 1, 2 * no_of_divisions + 1);

	for (int i = 1; i <= 2 * no_of_divisions + 1; i++)
		for (int j = 1; j <= 2 * no_of_divisions + 1; j++)
			transition_probability(i, j) = 0.0;

	// boundary values of the probabilities need to be entered
	transition_probability(1, 1) = 1.0 - uptick_prob;
	transition_probability(1, 2) = uptick_prob;
	transition_probability(2 * no_of_divisions + 1, 2 * no_of_divisions + 1) = 1.0 - downtick_prob;
	transition_probability(2 * no_of_divisions + 1, 2 * no_of_divisions) = downtick_prob;

	for (int i = 2; i <= 2 * no_of_divisions; i++)
		for (int j = 1; j <= 2 * no_of_divisions + 1; j++) 
		{
			if (j == (i - 1))
				transition_probability(i, j) = downtick_prob;
			if (j == i)
				transition_probability(i, j) = 1.0 - uptick_prob - downtick_prob;
			if (j == (i + 1))
				transition_probability(i, j) = uptick_prob;
		}

	Matrix V_t(2 * no_of_divisions + 1, 1);
	// value at expiration
	for (int i = 1; i <= 2 * no_of_divisions + 1; i++)
		V_t(i, 1) = max(0.0, strike_price - (initial_stock_price * pow(up_factor, i - 1 - no_of_divisions)));
	// value at intermediate stages
	Matrix discounted_one_step_forward_value(2 * no_of_divisions + 1, 1);
	Matrix value_if_option_is_exercised_now(2 * no_of_divisions + 1, 1);
	for (int i = no_of_divisions; i > 0; i--) 
	{
		// going backwards from expiration to zero-time
		discounted_one_step_forward_value = (transition_probability * V_t) / R;
		for (int j = 1; j <= 2 * no_of_divisions + 1; j++)
			value_if_option_is_exercised_now(j, 1) = max(0, strike_price - (initial_stock_price * pow(up_factor, j - 1 - no_of_divisions)));
		for (int j = 1; j <= 2 * no_of_divisions + 1; j++)
			V_t(j, 1) = max(value_if_option_is_exercised_now(j, 1), discounted_one_step_forward_value(j, 1));
	}
	return (V_t(no_of_divisions + 1, 1));
}

int main(int argc, char* argv[])
{
	sscanf_s(argv[1], "%lf", &expiration_time);
	sscanf_s(argv[2], "%d", &no_of_divisions);
	sscanf_s(argv[3], "%lf", &risk_free_rate);
	sscanf_s(argv[4], "%lf", &volatility);
	sscanf_s(argv[5], "%lf", &initial_stock_price);
	sscanf_s(argv[6], "%lf", &strike_price);

	up_factor = exp(volatility * sqrt(2 * expiration_time / ((float)no_of_divisions)));
	R = exp(risk_free_rate * expiration_time / ((float)no_of_divisions));
	uptick_prob = pow((sqrt(R) - (1 / sqrt(up_factor))) / (sqrt(up_factor) - (1 / sqrt(up_factor))), 2);
	downtick_prob = pow((sqrt(up_factor) - (sqrt(R))) / (sqrt(up_factor) - (1 / sqrt(up_factor))), 2);
	cout << "Recursive Trinomial American-Asian Option Pricing" << endl;
	cout << "Expiration Time (Years) = " << expiration_time << endl;
	cout << "Number of Divisions = " << no_of_divisions << endl;
	cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
	cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
	cout << "Initial Stock Price = " << initial_stock_price << endl;
	cout << "Strike Price = " << strike_price << endl;
	cout << "--------------------------------------" << endl;
	cout << "R = " << R << endl;
	cout << "Up Factor = " << up_factor << endl;
	cout << "Uptick Probability = " << uptick_prob << endl;
	cout << "Downtick Probability = " << downtick_prob << endl;
	cout << "Notick Probability = " << 1 - uptick_prob - downtick_prob << endl;
	cout << "--------------------------------------" << endl;
	clock_t start, end;
	start = clock();
	double call_price = american_call_option(0, 0, initial_stock_price);
	end = clock();
	cout << "Trinomial Price of an American Call Option (by Recursion) = " << call_price << endl;
	cout << "Running time = " << double(end - start) / CLOCKS_PER_SEC << "s" << std::endl;
	start = clock();
	double put_price = american_put_option(0, 0, initial_stock_price);
	end = clock();
	cout << "Trinomial Price of an American Put Option (by Recursion) = " << put_price << endl;
	cout << "Running time = " << double(end - start) / CLOCKS_PER_SEC << "s" << std::endl;
	start = clock();
	cout << "Trinomial Price of an American Call Option (by Dynamic Programming) = " << american_call_option_dyn_prog() << endl;
	end = clock();
	cout << "Running time = " << double(end - start) / CLOCKS_PER_SEC << "s" << std::endl;
	start = clock();
	cout << "Trinomial Price of an American Put Option (by Dynamic Programming) = " << american_put_option_dyn_prog() << endl;
	end = clock();
	cout << "Running time = " << double(end - start) / CLOCKS_PER_SEC << "s" << std::endl;
	cout << "--------------------------------------" << endl;
}