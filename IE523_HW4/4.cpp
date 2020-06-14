
// Computing the General Solution to Ax=b using the NEWMAT library

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>
#include <algorithm> 
#include "newmatap.h"           
#include "newmat.h"
#include "newmatio.h"
#include "minmax.h"
#define EPSILON 0.0000001
vector < vector < double > > solution_set;

int Rank(Matrix A)
{
	if (A.nrows() >= A.ncols())
	{
		DiagonalMatrix D(A.ncols());
		SVD(A, D);
		int rank = 0;
		for (int i = 1; i <= A.ncols(); i++)
			if (abs(D(i)) > EPSILON)
				rank++;
		return (rank);
	}
	else {
		DiagonalMatrix D(A.nrows());
		SVD(A.t(), D);
		int rank = 0;
		for (int i = 1; i <= A.nrows(); i++)
			if (abs(D(i)) > EPSILON)
				rank++;
		return (rank);
	}
}

// This routine finds the basic solution where we have rank(A) many columns and y is in the linear span of the cols of A
ColumnVector Find_Basic_Solution(Matrix A, ColumnVector y)
{
	//A should be of full rank and y is a column vector
	if (Rank(A) == min(A.ncols(), A.nrows()))
	{
		return ((A.t() * A).i()) * A.t() * y;
	}
	else {
		cout << "Matrix: " << endl;
		cout << setw(9) << setprecision(3) << A;
		cout << "is not of full-rank... exiting..." << endl;
		exit(0);
	}
}

void create_solution_vector(Matrix A, ColumnVector y)
{    
	int rank_of_A = Rank(A);

	{
		vector <bool> did_i_pick_this_col(A.Ncols());
		fill(did_i_pick_this_col.begin(), 
			 did_i_pick_this_col.begin() + rank_of_A, true);
		do 
		{
            vector <vector <double> > B;
			int k = 0;
			for (int i = 0; i < A.ncols(); i++) 
			{
				// fill this code appropriately
				if (did_i_pick_this_col[i] == true)
				{
					B.push_back(vector <double>());//col
					for (int j = 1; j <= A.nrows(); j++) 
					{
						B[k].push_back(A(j, i + 1));
					}
					k++;
				}
			}
			Matrix Column_reduced_version_of_A(A.nrows(), rank_of_A);
			for (int j = 1; j <= A.nrows(); j++)
				for (int k = 1; k <= rank_of_A; k++)
					Column_reduced_version_of_A(j,k) = B[k-1][j-1];
			if (Rank(Column_reduced_version_of_A) == rank_of_A)
			{
				ColumnVector soln(rank_of_A);
				soln = Find_Basic_Solution(Column_reduced_version_of_A, y);
				int k = 1;
                vector <double> temp;
				for (int i = 0; i < A.ncols(); ++i)
				{
					if (!did_i_pick_this_col[i])
					{
						temp.push_back(0);
					}
					else
					{
						temp.push_back(soln(k));
						k++;
					}
				}
                solution_set.push_back(temp);
            }
		} while (prev_permutation(did_i_pick_this_col.begin(), did_i_pick_this_col.end())); 
	}
	
}

void trim_the_solution_set(Matrix A, ColumnVector y)
{
	Matrix solution(A.ncols(), solution_set.size());
	for (int j = 1; j <= A.ncols(); j++)
		for (int k = 1; k <= solution_set.size(); k++)
			solution(j,k) = solution_set[k-1][j-1];
	int rank_of_solution_set = Rank(solution);
	
	{
		vector <bool> did_i_pick_this(solution.ncols());
		fill(did_i_pick_this.begin(), 
			 did_i_pick_this.begin() + rank_of_solution_set, true);
		do 
		{
			vector <vector <double> > B;
			int k = 0;
			for (int i = 0; i < solution.ncols(); i++)
			{
				if (did_i_pick_this[i])
				{
					B.push_back(vector <double>());
					for (int j = 1; j <= solution.nrows(); j++)
					{
						B[k].push_back(solution(j, i + 1));
					}
					k++;
				}
			}
			Matrix Column_reduced_version(solution.nrows(), rank_of_solution_set);
			for (int j = 1; j <= solution.nrows(); j++)
				for (int k = 1; k <= rank_of_solution_set; k++)
					Column_reduced_version(j,k) = B[k-1][j-1];
			if (Rank(Column_reduced_version) == rank_of_solution_set)
			{
				cout << "Solving:" << endl;
				cout << setw(10) << A << endl;
				cout << "* solution = " << endl << endl;
				cout << setw(10) << y << endl;
				cout << "solution is the affine-combination of these vectors" << endl;
				cout << setw(10) <<  Column_reduced_version << endl;
				cout << "Verification: (check each column below and y-vector above)" << endl;
				cout << setw(10) << A*Column_reduced_version << endl;
				exit(0);
			}
		} while (prev_permutation(did_i_pick_this.begin(), did_i_pick_this.end())); 
	}
}

void process_data(char* argv[])
{
	int no_of_rows, no_of_cols;
	double element;
	ifstream input_file(argv[1]);
	if (input_file.is_open())
	{
		input_file >> no_of_rows; 
		input_file >> no_of_cols; 
		
		// Create the space for the matrix and column vector
		Matrix A(no_of_rows, no_of_cols);
		ColumnVector y(no_of_rows);
		
		//read the A matrix
		for (int i = 1; i <= no_of_rows; i++)
		{
			for (int j = 1; j <= no_of_cols; j++)
			{
				input_file >> element;
				A(i,j) = element;
			}
		}
		
		// read the y vector
		for (int i = 1; i <= no_of_rows; i++)
		{
			input_file >> element;
			y(i) = element;
		}
		
		// check if there is a solution, if "yes" then compute genl soln
		Matrix B(no_of_rows, no_of_cols+1);
		B = (A | y);
		if (Rank(B) == Rank(A))
		{
			create_solution_vector(A, y);
			trim_the_solution_set(A, y);
		}
		else 
		{
			// we found a submatrix of solution with the same rank
			cout << "Solving:" << endl;
			cout << setw(10) << A << endl;
			cout << "* vector_x = " << endl << endl;
			cout << setw(10) << y << endl;
			cout << "There is no solution to this equation" << endl;
			exit(0);
		}
	}
	else 
	{
		cout << "Input file: " << argv[1] << " does not exist in present folder"  << endl;
	}	
}
	
int main (int argc, char* argv[])
{
	if (argc == 1) 
		cout << "Input filename missing" << endl;
	else 
		process_data(argv);
}


	
	