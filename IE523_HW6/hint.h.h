
#ifndef ALICE_AND_BOB
#define ALICE_AND_BOB
#include <cmath>
#include <vector>
#include <fstream>
#include <cstdlib>
using namespace std;

class I_have_nothing_apropos_for_this_class
{
private:
	double alice_probability, bob_probability;
	double get_uniform()
	{
		return (((double)rand()) / (pow(2.0, 15.0) - 1.0));
	}

	int take(int n, int i)
	{
		if (i != 0) 
		{
			return ((take(n, i - 1) * (n - i + 1)) / i);
		}
		else 
		{
			return 1;
		}
	}
	
	double theoretical_value(double q, double p, int n)
	{
		double f,w;
		f = 0;
		for (int r = 0; r <= n; r++)
		{
			w = 0;
			for (int s = r + 1; s <= n; s++)
			{
				w = w + take(n, s) * pow(q, s) * pow(1 - q, n - s) ;
			}
			f = f + take(n, r) * pow(p, r) * pow(1 - p, n - r) * w;
		}
		return f;
	}

public: 
	void set_probability(double alice_p, double bob_p)
	{
		alice_probability = alice_p;
		bob_probability = bob_p;
	}
	
	// probability of Alice winning the game.
	double simulated_value(int number_of_coin_tosses_in_each_game, int no_of_trials)
	{
		int no_of_wins_for_alice = 0;
		for (int i = 0; i < no_of_trials; i++) 
		{
			int number_of_heads_for_alice = 0;
			int number_of_heads_for_bob = 0;
			for (int j = 0; j < number_of_coin_tosses_in_each_game; j++) 
			{
				if (get_uniform() < alice_probability) 
					number_of_heads_for_alice++;
				if (get_uniform() < bob_probability)
					number_of_heads_for_bob++;
			}
			if (number_of_heads_for_alice > number_of_heads_for_bob)
				no_of_wins_for_alice++;
		}
		return (((double) no_of_wins_for_alice)/((double) no_of_trials));
	}
		
	int search_result()
	{
		int n;
		double d = 0;
		for (n = 1; n <= 100; n++)
		{
			d = theoretical_value(alice_probability, bob_probability, n+1) - theoretical_value(alice_probability, bob_probability, n);
			if (d <= 0)
				return n;
		}
	}

	vector <double> simulation;
	vector <double> theory;
	void output() 
	{
		ofstream file1;
		ofstream file2;
		file1.open("data_of_simulation");
		file2.open("data_of_theory");
		for (int i = 0; i < 30; i++) 
		{
			theory.push_back(theoretical_value(alice_probability, bob_probability, i + 1));
			simulation.push_back(simulated_value(i + 1, 500000));
			file1 << simulation[i] << endl;
			file2 << theory[i] << endl;
		}
	}
};
#endif









