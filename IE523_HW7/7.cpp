
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>
#include <algorithm> 
#include <time.h>
#include "newmatap.h"           
#include "newmat.h"
#include "newmatio.h"
#include "minmax.h"

Matrix repeated_squaring(Matrix A, int exponent, int dimension)
{
	if (exponent <= 0)
		return IdentityMatrix(dimension);
	else if (exponent % 2 == 1)
	{
		return A * repeated_squaring(A * A, (exponent - 1) / 2, dimension);
	}
	else
	{
		return repeated_squaring(A * A, exponent / 2, dimension);
	}
}

Matrix direct_multiplying(Matrix A, int exponent)
{
	Matrix C;
	C = A;
	for (int i = 1; i <= exponent - 1; i++)
	{
		C = A * C;
	}
	return C;
}

void output(int dimension)
{
	ofstream file;
	file.open("time");
	clock_t start, end;
	float time_R, time_D;
	Matrix A(dimension, dimension);
	for (int i = 1; i <= dimension; i++)
	{
		for (int j = 1; j <= dimension; j++)
		{
			A(i, j) = (double)rand() / (pow(2.0, 15.0) - 1.0) * 10 - 5;
		}
	}
	for (int j = 50; j <= 1000; j += 50)
	{

		start = clock();
		repeated_squaring(A, j, dimension);
		end = clock();
		time_R = ((float)end - (float)start) / CLOCKS_PER_SEC;

		start = clock();
		direct_multiplying(A, j);
		end = clock();
		time_D = ((float)end - (float)start) / CLOCKS_PER_SEC;

		file << j << "," << time_R << "," << time_D << endl;
	}
}

void process(int n,int k)
{
	// Create the space for the matrix
	Matrix A(n, n);
	for (int i = 1; i <= n; i++)
	{
		for (int j = 1; j <= n; j++)
		{
			A(i, j) = (double)rand() / (pow(2.0, 15.0) - 1.0) * 10 - 5;
		}
	}
	cout << "The number of rows/columns in the square matrix is: " << n << endl;
	cout << "The exponent is: " << k << endl;

	clock_t time_before1, time_after1, time_before2, time_after2;

	// Computing by Repeated Squaring Algorithm
	Matrix B(n, n);
	time_before1 = clock();
	B = repeated_squaring(A, k, n);
	time_after1 = clock();
	float diff1 = ((float)time_after1 - (float)time_before1);
	cout << "Repeated Squaring Result:" << endl;
	cout << "It took " << diff1 / CLOCKS_PER_SEC << " seconds to complete" << endl;

	// Computing by Direct Multiplication
	Matrix C(n, n);
	time_before2 = clock();
	C = direct_multiplying(A, k);
	time_after2 = clock();
	float diff2 = ((float)time_after2 - (float)time_before2);
	cout << "Direct Multiplication Result:" << endl;;
	cout << "It took " << diff2 / CLOCKS_PER_SEC << " seconds to complete" << endl;

	exit(0);
}

int main(int argc, char* argv[])
{
	int n, k;
	sscanf_s(argv[1], "%d", &k);
	sscanf_s(argv[2], "%d", &n);
	process(n, k);
}

