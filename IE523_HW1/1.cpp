// N Queens Problem via (Backtracking, which is implemented by) Recursion 
#include <iostream>
#include "Part 1.h.h"
#include "Part 2.h.h"

int main (int argc, char ** argv) 
{
    Board x;
    int board_size;
	cin>> board_size;
	x.nQueens(board_size);
    return 0;
}
