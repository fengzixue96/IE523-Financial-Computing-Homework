
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <time.h>
#include "hint.h.h"
using namespace std;
	
int main (int argc, char* argv[])
{
	I_have_nothing_apropos_for_this_class x;
	double alice_success_prob, bob_success_prob;
	sscanf_s(argv[1], "%lf", &alice_success_prob);
	sscanf_s(argv[2], "%lf", &bob_success_prob);
	//cin >> alice_success_prob >> bob_success_prob;
	
	cout << "Probability of success for Alice = " << alice_success_prob << endl;
	cout << "Probability of success for Bob = " << bob_success_prob << endl;
	
	x.set_probability(alice_success_prob, bob_success_prob);

	clock_t start, end;
	start = clock();
	int optimal = x.search_result();
	end = clock();
	if (optimal > 0)
	{
		cout << "The optimal number of coin tosses in each game is " << optimal << endl;
		cout << "This program took me " << double(end - start) / CLOCKS_PER_SEC << " seconds" << endl;
	}
	else {
		cout << "The optimal number of coin tosses in each game exceeds 100... Quitting" << endl;
	}
	x.output();
}
	
	
	
