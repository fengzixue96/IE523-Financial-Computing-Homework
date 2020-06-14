
// Soduku Solver using Brute-Force Search implemted using recursion.

#include <iostream>
#include "Assignment2.h.h"

int main (int argc, char * const argv[]) 
{
	Sudoku x;
	x.read_puzzle(argc, argv);
	x.print_puzzle();
	x.Solve(0,0);

    return 0;
}
