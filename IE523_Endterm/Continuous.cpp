
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <random>
#include <chrono>
#include "normdist.h"
using namespace std;

double risk_free_rate, strike_price, initial_stock_price, expiration_time, volatility, barrier_price;
int no_of_trials, no_of_divisions;

unsigned seed = (unsigned)std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);

double get_uniform()
{
	std::uniform_real_distribution <double> distribution(0, 1.0);
	double number = distribution(generator);
	return (number);
}

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

double option_price_call_black_scholes(const double& S,       // spot (underlying) price
	const double& K,       // strike (exercise) price,
	const double& r,       // interest rate
	const double& sigma,   // volatility
	const double& time)      // time to maturity
{
	double time_sqrt = sqrt(time);
	double d1 = (log(S / K) + r * time) / (sigma * time_sqrt) + 0.5 * sigma * time_sqrt;
	double d2 = d1 - (sigma * time_sqrt);
	return S * N(d1) - K * exp(-r * time) * N(d2);
};

double option_price_put_black_scholes(const double& S,      // spot price
	const double& K,      // Strike (exercise) price,
	const double& r,      // interest rate
	const double& sigma,  // volatility
	const double& time)
{
	double time_sqrt = sqrt(time);
	double d1 = (log(S / K) + r * time) / (sigma * time_sqrt) + 0.5 * sigma * time_sqrt;
	double d2 = d1 - (sigma * time_sqrt);
	return K * exp(-r * time) * N(-d2) - S * N(-d1);
};

double closed_form_down_and_out_european_call_option()
{
	// I took this formula from Wilmott, Howison and Dewynne, "The Mathematics of Financial Derivatives"
	double K = (2 * risk_free_rate) / (volatility * volatility);
	double A = option_price_call_black_scholes(initial_stock_price, strike_price,
		risk_free_rate, volatility, expiration_time);
	double B = (barrier_price * barrier_price) / initial_stock_price;
	double C = pow(initial_stock_price / barrier_price, -(K - 1));
	double D = option_price_call_black_scholes(B, strike_price, risk_free_rate, volatility, expiration_time);
	return (A - D * C);
}

double closed_form_down_and_in_european_put_option()
{
	// just making it easier by renaming the global variables locally
	double S = initial_stock_price;
	double r = risk_free_rate;
	double T = expiration_time;
	double sigma = volatility;
	double H = barrier_price;
	double X = strike_price;

	// Took these formulae from some online reference
	double lambda = (r + ((sigma * sigma) / 2)) / (sigma * sigma);
	double temp = 2 * lambda - 2.0;
	double x1 = (log(S / H) / (sigma * sqrt(T))) + (lambda * sigma * sqrt(T));
	double y = (log(H * H / (S * X)) / (sigma * sqrt(T))) + (lambda * sigma * sqrt(T));
	double y1 = (log(H / S) / (sigma * sqrt(T))) + (lambda * sigma * sqrt(T));
	return (-S * N(-x1) + X * exp(-r * T) * N(-x1 + sigma * sqrt(T)) +
		S * pow(H / S, 2 * lambda) * (N(y) - N(y1)) -
		X * exp(-r * T) * pow(H / S, temp) * (N(y - sigma * sqrt(T)) - N(y1 - sigma * sqrt(T))));
}

double closed_form_down_and_out_european_put_option()
{
	double vanilla_put = option_price_put_black_scholes(initial_stock_price, strike_price,
		risk_free_rate, volatility, expiration_time);
	double put_down_in = closed_form_down_and_in_european_put_option();
	return (vanilla_put - put_down_in);
}

int main(int argc, const char* argv[]) {
	sscanf_s(argv[1], "%lf", &expiration_time);
	sscanf_s(argv[2], "%lf", &risk_free_rate);
	sscanf_s(argv[3], "%lf", &volatility);
	sscanf_s(argv[4], "%lf", &initial_stock_price);
	sscanf_s(argv[5], "%lf", &strike_price);
	sscanf_s(argv[6], "%d", &no_of_trials);
	sscanf_s(argv[7], "%d", &no_of_divisions);
	sscanf_s(argv[8], "%lf", &barrier_price);

	const double delta_T = expiration_time / ((double)no_of_divisions);
	const double delta_R = (risk_free_rate - 0.5 * pow(volatility, 2)) * delta_T;
	const double delta_SD = volatility * sqrt(delta_T);

	cout << "European Down-and-out Continuous Barrier Options Pricing via Monte Carlo Simulation" << endl;
	cout << "--------------------------------" << endl;
	cout << "Expiration Time (Years) = " << expiration_time << endl;
	cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
	cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
	cout << "Initial Stock Price = " << initial_stock_price << endl;
	cout << "Strike Price = " << strike_price << endl;
	cout << "Barrier Price = " << barrier_price << endl;
	cout << "Number of Trials = " << no_of_trials << endl;
	cout << "Number of Divisions = " << no_of_divisions << endl;
	cout << "--------------------------------" << endl;

	double put_option_price = 0.0;
	double call_option_price = 0.0;
	double put_option_price_adj = 0.0;
	double call_option_price_adj = 0.0;
	for (int i = 0; i < no_of_trials; i++)
	{
		double current_stock_price1 = initial_stock_price;
		double current_stock_price2 = initial_stock_price;
		double current_stock_price3 = initial_stock_price;
		double current_stock_price4 = initial_stock_price;

		double current_stock_price1_adj = initial_stock_price;
		double current_stock_price2_adj = initial_stock_price;
		double current_stock_price3_adj = initial_stock_price;
		double current_stock_price4_adj = initial_stock_price;

		for (int j = 0; j < no_of_divisions; j++)
		{
			// create the unit normal variates using the Box-Muller Transform
			double x = get_uniform();
			double y = get_uniform();
			double a = sqrt(-2.0 * log(x)) * cos(6.283185307999998 * y);
			double b = sqrt(-2.0 * log(x)) * sin(6.283185307999998 * y);

			// check if the current_stock_price > barrier price ?

			current_stock_price1 = ((current_stock_price1 * exp(delta_R + delta_SD * a)) > barrier_price) ? (current_stock_price1 * exp(delta_R + delta_SD * a)) : 0;
			current_stock_price2 = ((current_stock_price2 * exp(delta_R - delta_SD * a)) > barrier_price) ? (current_stock_price2 * exp(delta_R - delta_SD * a)) : 0;
			current_stock_price3 = ((current_stock_price3 * exp(delta_R + delta_SD * b)) > barrier_price) ? (current_stock_price3 * exp(delta_R + delta_SD * b)) : 0;
			current_stock_price4 = ((current_stock_price4 * exp(delta_R - delta_SD * b)) > barrier_price) ? (current_stock_price4 * exp(delta_R - delta_SD * b)) : 0;

			current_stock_price1_adj = current_stock_price1_adj * exp(delta_R + delta_SD * a);
			current_stock_price2_adj = current_stock_price2_adj * exp(delta_R - delta_SD * a);
			current_stock_price3_adj = current_stock_price3_adj * exp(delta_R + delta_SD * b);
			current_stock_price4_adj = current_stock_price4_adj * exp(delta_R - delta_SD * b);
		}
		double pd1 = (current_stock_price1 != 0) ? 1 : 0;
		double pd2 = (current_stock_price2 != 0) ? 1 : 0;
		double pd3 = (current_stock_price3 != 0) ? 1 : 0;
		double pd4 = (current_stock_price4 != 0) ? 1 : 0;

		call_option_price += (max(0.0, current_stock_price1 - strike_price) * pd1 + max(0.0, current_stock_price2 - strike_price) * pd2 + max(0.0, current_stock_price3 - strike_price) * pd3 + max(0.0, current_stock_price4 - strike_price) * pd4) / 4.0;
		put_option_price += (max(0.0, strike_price - current_stock_price1) * pd1 + max(0.0, strike_price - current_stock_price2) * pd2 + max(0.0, strike_price - current_stock_price3) * pd3 + max(0.0, strike_price - current_stock_price4) * pd4) / 4.0;

		double pc1, pc2, pc3, pc4;
		pc1 = (initial_stock_price <= barrier_price || current_stock_price1_adj <= barrier_price) ? 1.0 : (exp((-2.0 / (volatility * volatility * expiration_time)) * log(initial_stock_price / barrier_price) * log(current_stock_price1_adj / barrier_price)));
		pc2 = (initial_stock_price <= barrier_price || current_stock_price2_adj <= barrier_price) ? 1.0 : (exp((-2.0 / (volatility * volatility * expiration_time)) * log(initial_stock_price / barrier_price) * log(current_stock_price2_adj / barrier_price)));
		pc3 = (initial_stock_price <= barrier_price || current_stock_price3_adj <= barrier_price) ? 1.0 : (exp((-2.0 / (volatility * volatility * expiration_time)) * log(initial_stock_price / barrier_price) * log(current_stock_price3_adj / barrier_price)));
		pc4 = (initial_stock_price <= barrier_price || current_stock_price4_adj <= barrier_price) ? 1.0 : (exp((-2.0 / (volatility * volatility * expiration_time)) * log(initial_stock_price / barrier_price) * log(current_stock_price4_adj / barrier_price)));

		call_option_price_adj += (max(0.0, current_stock_price1_adj - strike_price) * (1.0 - pc1) + max(0.0, current_stock_price2_adj - strike_price) * (1.0 - pc2) + max(0.0, current_stock_price3_adj - strike_price) * (1.0 - pc3) + max(0.0, current_stock_price4_adj - strike_price) * (1.0 - pc4)) / 4.0; 
		put_option_price_adj += (max(0.0, strike_price - current_stock_price1_adj) * (1.0 - pc1) + max(0.0, strike_price - current_stock_price2_adj) * (1.0 - pc2) + max(0.0, strike_price - current_stock_price3_adj) * (1.0 - pc3) + max(0.0, strike_price - current_stock_price4_adj) * (1.0 - pc4)) / 4.0;

	}
	call_option_price = exp(-risk_free_rate * expiration_time) * (call_option_price / ((double)no_of_trials));
	put_option_price = exp(-risk_free_rate * expiration_time) * (put_option_price / ((double)no_of_trials));
	call_option_price_adj = exp(-risk_free_rate * expiration_time) * (call_option_price_adj / ((double)no_of_trials));
	put_option_price_adj = exp(-risk_free_rate * expiration_time) * (put_option_price_adj / ((double)no_of_trials));

	cout << "--------------------------------" << endl;
	cout << "The average Call Price by explict simulation = " << call_option_price << endl;
	cout << "The Call Price using (1-p)-adjustment term   = " << call_option_price_adj << endl;
	cout << "Theoretical Call Price                       = " << closed_form_down_and_out_european_call_option() << endl;

	cout << "--------------------------------" << endl;
	cout << "The average Put Price by explict simulation = " << put_option_price << endl;
	cout << "The Put Price using (1-p)-adjustment term   = " << put_option_price_adj << endl;
	cout << "Theoretical Call Price                      = " << closed_form_down_and_out_european_put_option() << endl;

	return 0;
}