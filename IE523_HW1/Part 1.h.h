
// N Queens Problem
// Just Print One Solution

#ifndef N_queens
#define N_queens
#include <iostream>
using namespace std;
class Board
{
    
    // private data member: size of the board
    int size;
    
    // pointer-to-pointer initialization of the board
    int **chess_board;
    
    // private member function:  returns 'false' if the (row, col) position is not safe.
    bool is_this_position_safe(int row, int col)
    {
		bool result = true;
		int i, j,k;
		//check column
		for (i = 0; i < size; i++)
		{
			if (chess_board[row][i] == 1)
			{
				result = false;
			}
		}
		//check top left corner 
		for (j = 0; j<=row && j<=col; j++)
		{
			if (chess_board[row-j][col-j]==1)
			{
				return false;
			}
		}
		//check left bottom
		for (k = 0; k < size-row && k <= col; k++)
		{
			if (chess_board[row + k][col - k]==1)
			{
				return false;
			}
		}
		return (result);
        // returns "true" if the (row,col) position is safe.  
		//If it is unsafe(i.e. some other queen can threaten this position), return "false"
    }
    
    // private member function: initializes the (n x n) chessboard
    void initialize(int n)
    {
        size = n;
		// chess_board is a pointer-to-pointer
		chess_board = new int* [size];
		for (int i = 0; i < size; i++)
		{
			chess_board[i] = new int[size];
		}
		// method to initialize the (n x n) chessboard.  
		for (int j=0; j < size; j++)
		{
			for (int k=0; k < size; k++)
			{
				chess_board[j][k] = 0;
			}
		}
		//Once initialized, put zeros in all entries.  
		//Later on, if you placed a queen in the (i,j)-th position, then chessboard[i][j] will be 1.
    }
    
    // private member function: prints the board position
    void print_board()
    {
        cout << "\n" << size << "-Queens Problem Solution" <<  endl;
		for (int k = 0; k < 2*size+1; k++)
		{
			cout << "-";
		}
		cout <<  endl;
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				if (chess_board[i][j] == 0)
					cout << " -";
				else
					cout << " Q";
			}
			cout << endl;
		}
		for (int k = 0; k < 2 * size + 1; k++)
		{
			cout << "-";
		}
		cout << endl;
    }
    
    // private member function: recursive backtracking
	// pseudocode format in figure 1 of the description of the first
    bool solve(int col)
	{
		if (col == size)
			return true;
		else
			{
				for (int i = 0; i < size; i++)
				{
					if (is_this_position_safe(i, col))
					{
						chess_board[i][col] = 1;
						if (solve(col + 1))
							return true;
						else
							chess_board[i][col] = 0;
					}
				}
			}
        return false;
	}
    
public:
    // Solves the n-Queens problem by (recursive) backtracking
    void nQueens(int n)
    {
		initialize(n);    
		if(solve(0))
			print_board();
		else
			cout << "There is no solution to the " << n << "-Queens Problem." << endl;
    }
};
#endif
