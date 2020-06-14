
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <time.h>
using namespace std;

double up_factor, uptick_prob, downtick_prob, risk_free_rate, strike_price;
double initial_stock_price, expiration_time, volatility, R;
int no_of_divisions;
float** dic;

//Pricing an Option using Black-Scholes
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

double option_price_put_black_scholes(const double& S,      // spot price
	const double& K,  
	const double& r,  
	const double& sigma, // volatility 
	const double& time)  // time to maturity 
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

//Pricing an American/European-Option using a (memorized) Trinomial Model by Recursion
void initialize()
{
	dic = new float* [no_of_divisions+1];
	for (int i = 0; i < no_of_divisions+1; i++)
	{
		dic[i] = new float[no_of_divisions*2+1];
	}
	for (int j = 0; j < no_of_divisions+1; j++)
	{
		for (int k = 0; k < no_of_divisions * 2 + 1; k++)
		{
			dic[j][k] = -1;
		}
	}
}

double european_call_option(int k, int i)
{
	if (dic[k][i + no_of_divisions] != -1)
	{
		return dic[k][i + no_of_divisions];
	}
	else 
	{
		if (k == no_of_divisions)
		{
			dic[k][i + no_of_divisions] = max(0.0, (initial_stock_price * pow(up_factor, ((double)i))) - strike_price);
			return dic[k][i + no_of_divisions];
		}
		else
		{
			dic[k][i + no_of_divisions] = (uptick_prob * european_call_option(k + 1, i + 1) +
				(1 - uptick_prob - downtick_prob) * european_call_option(k + 1, i) +
				downtick_prob * european_call_option(k + 1, i - 1)) / R;
			return dic[k][i + no_of_divisions];
		}
	}
}

double european_put_option(int k, int i) 
{
	if (dic[k][i + no_of_divisions] != -1)
	{
		return dic[k][i + no_of_divisions];
	}
	else
	{
		if (k == no_of_divisions)
		{
			dic[k][i + no_of_divisions] = max(0.0, strike_price - (initial_stock_price * pow(up_factor, ((double)i))));
			return dic[k][i + no_of_divisions];
		}
		else
		{
			dic[k][i + no_of_divisions] = (uptick_prob * european_put_option(k + 1, i + 1) +
				(1 - uptick_prob - downtick_prob) * european_put_option(k + 1, i) +
				downtick_prob * european_put_option(k + 1, i - 1)) / R;
			return dic[k][i + no_of_divisions];
		}
	}
}

double american_call_option(int k, int i, double current_stock_price)
{
	if(dic[k][i + no_of_divisions]!=-1)
	{
		return dic[k][i + no_of_divisions];
	}
	else 
		if(k == no_of_divisions)
		{
			dic[k][i + no_of_divisions] = max(0.0, (current_stock_price - strike_price));
			return dic[k][i + no_of_divisions];
		}
		else
		{
			dic[k][i + no_of_divisions] = max((current_stock_price - strike_price),
				(uptick_prob * american_call_option(k + 1, i + 1, current_stock_price * up_factor) +
				(1 - uptick_prob - downtick_prob) * american_call_option(k + 1, i, current_stock_price) +
				downtick_prob * american_call_option(k + 1, i - 1, current_stock_price / up_factor)) / R);
			return dic[k][i + no_of_divisions];
		}
}

double american_put_option(int k, int i, double current_stock_price)
{
	if (dic[k][i + no_of_divisions] != -1)
	{
		return dic[k][i + no_of_divisions];
	}
	else
	{
		if (k == no_of_divisions)
		{
			dic[k][i + no_of_divisions] = max(0.0, (strike_price - current_stock_price));
			return dic[k][i + no_of_divisions];
		}
		else
		{
			dic[k][i + no_of_divisions] = max((strike_price - current_stock_price),
				(uptick_prob * american_put_option(k + 1, i + 1, current_stock_price * up_factor) +
				(1 - uptick_prob - downtick_prob) * american_put_option(k + 1, i, current_stock_price) +
					downtick_prob * american_put_option(k + 1, i - 1, current_stock_price / up_factor)) / R);
			return dic[k][i + no_of_divisions];
		}
	}

}

int main(int argc, char* argv[])
{
	sscanf_s(argv[1], "%lf", &expiration_time);
	sscanf_s(argv[2], "%d", &no_of_divisions);
	sscanf_s(argv[3], "%lf", &risk_free_rate);
	sscanf_s(argv[4], "%lf", &volatility);
	sscanf_s(argv[5], "%lf", &initial_stock_price);
	sscanf_s(argv[6], "%lf", &strike_price);
	clock_t start, end;
	
	up_factor = exp(volatility * sqrt(2 * expiration_time / ((float)no_of_divisions)));
	R = exp(risk_free_rate * expiration_time / ((float)no_of_divisions));
	uptick_prob = pow((sqrt(R) - (1 / sqrt(up_factor))) / (sqrt(up_factor) - (1 / sqrt(up_factor))), 2);
	downtick_prob = pow((sqrt(up_factor) - (sqrt(R))) / (sqrt(up_factor) - (1 / sqrt(up_factor))), 2);
	
	cout << "(Memoized) Recursive Trinomial European Option Pricing" << endl;
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

	initialize();
	start = clock();
	double call_price_e = european_call_option(0, 0);
	end = clock();
	cout << "Trinomial Price of an European Call Option = " << call_price_e << endl;
	cout << "Call Price according to Black-Scholes = " <<
			option_price_call_black_scholes(initial_stock_price, strike_price, risk_free_rate,
			volatility, expiration_time) << endl;
	cout << "Running time = " << double(end - start) / CLOCKS_PER_SEC << "s" << std::endl;
	cout << "--------------------------------------" << endl;
	initialize();
	start = clock();
	double put_price_e = european_put_option(0, 0);
	end = clock();
	cout << "Trinomial Price of an European Put Option = " << put_price_e << endl;
	cout << "Put Price according to Black-Scholes = " <<
			option_price_put_black_scholes(initial_stock_price, strike_price, risk_free_rate,
			volatility, expiration_time) << endl;
	cout << "Running time = " << double(end - start) / CLOCKS_PER_SEC << "s" << std::endl;
	cout << "--------------------------------------" << endl;
	cout << "Verifying Put-Call Parity: S+P-C = Kexp(-r*T)" << endl;
	cout << initial_stock_price << " + " << put_price_e << " - " << call_price_e;
	cout << " = " << strike_price << "exp(-" << risk_free_rate << " * " << expiration_time << ")" << endl;
	cout << initial_stock_price + put_price_e - call_price_e << " = " << strike_price * exp(-risk_free_rate * expiration_time) << endl;
	cout << "--------------------------------------" << endl;

	cout << "(Memoized) Recursive Trinomial American Option Pricing" << endl;
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

	initialize();
	start = clock();
	double call_price_a = american_call_option(0, 0, initial_stock_price);
	end = clock();
	cout << "Trinomial Price of an American Call Option = " << call_price_a << endl;
	cout << "Running time = " << double(end - start) / CLOCKS_PER_SEC << "s" << std::endl;
	initialize();
	start = clock();
	double put_price_a = american_put_option(0, 0, initial_stock_price);
	end = clock();
	cout << "Trinomial Price of an American Put Option = " << put_price_a << endl;
	cout << "Running time = " << double(end - start) / CLOCKS_PER_SEC << "s" << std::endl;
	cout << "--------------------------------------" << endl;
	cout << "Let us verify the Put-Call Parity: S+P-C = Kexp(-r*T) for American Options" << endl;
	cout << initial_stock_price << " + " << put_price_a << " - " << call_price_a;
	cout << " = " << strike_price << "exp(-" << risk_free_rate << " * " << expiration_time << ")" << endl;
	cout << initial_stock_price + put_price_a - call_price_a << " ?=? " << strike_price * exp(-risk_free_rate * expiration_time) << endl;
	if (abs(initial_stock_price + put_price_a - call_price_a - strike_price * exp(-risk_free_rate * expiration_time)) <= 1e-3)
		cout << "Looks like Put-Call Parity holds within three decimal places" << endl;
	else
		cout << "Looks like Put-Call Parity does NOT hold" << endl;
	cout << "--------------------------------------" << endl;
}