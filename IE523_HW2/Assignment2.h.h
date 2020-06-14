#ifndef sudoku
#define sudoku

#include <iostream>
#include <vector>
#include <fstream>
using std::vector;
using namespace std;
class Sudoku
{
	// Private
	int puzzle[9][9];
	int sum = 0;

	// Private member function that checks if the named row is valid
	bool row_valid(int row,int a)
	{
		for (int i = 0; i < 9; i++)
			if (puzzle[row][i]==a)
				return false;
		return true;
	}

	// Private member function that checks if the named column is valid
	bool col_valid(int col,int a)
	{
		for (int i = 0; i < 9; i++)
			if (puzzle[i][col]==a)
				return false;
		return true;
	}

	// Private member function that checks if the named 3x3 block is valid
	bool block_valid(int row, int col,int a)
	{
		int i = row - (row % 3);
        int j = col - (col % 3);
		for (int k = 0; k < 3; k++)
			for (int h = 0; h < 3; h++)
				if (puzzle[i + k][j + h] == a)
					return false;
		return true;
	}

public:
	void read_puzzle(int argc, char* const argv[])
	{
		ifstream input_file(argv[1]);
		if (input_file.is_open())
		{
			for (int i = 0; i < 9; i++)
				for (int j = 0; j < 9; j++)
					input_file >> puzzle[i][j];
			input_file.close();
		}
		else
			cout << "The file is not exist." << endl;
	}

	void print_puzzle()
	{
		if (sum == 0)
			cout << endl << "Board Position" << endl;
		else
			cout << endl << "Solution #"<<sum<<endl<<"Board Position" << endl;
		for (int i = 0; i < 9; i++)
		{
			for (int j = 0; j < 9; j++)
			{
				// check if we have a legitimate integer between 1 and 9
				if ((puzzle[i][j] >= 1) && (puzzle[i][j] <= 9))
				{
					// printing initial value of the puzzle with some formatting
					cout << puzzle[i][j] << " ";
				}
				else {
					// printing initial value of the puzzle with some formatting
					cout << "X ";
				}
			}
			cout << endl;
		}
	}

	// Public member function that (recursively) implements the brute-force 
	// search for possible solutions to the incomplete Sudoku puzzle
	bool Solve(int row, int col)
	{
		int i, j;
		for (int i = 0; i < 9; i++)
		{
			for (int j = 0; j < 9; j++)
			{
				if (puzzle[i][j] == 0)
				{
					for (int k = 1; k < 10; k++)
					{
						if (block_valid(i, j, k) && col_valid(j,k) && row_valid(i,k))
						{
							puzzle[i][j] = k;
							if (Solve(i, j))
								return true;
							else
								puzzle[i][j] = 0;
						}
					}
					return false;
				}
			}
		}
		sum = sum + 1;
		print_puzzle();
		return false;
	}
};

#endif