
#include <iostream>
#include <time.h>
#include <stdio.h>
using namespace std;
#include <algorithm>
class Board
{
	float num = 0;
	float** dic;
	float seek(int a, int b)
	{
		if (dic[a][b] != -1)
		{
			return dic[a][b];
		}
		else
		{
			return -1;
		}
	}
	void in(int a,int b,float c)
	{
		dic[a][b] = c;
	}
	void initialize()
	{
		dic = new float* [100];
		for (int i = 0; i < 100; i++)
		{
			dic[i] = new float[100];
		}
		for (int j = 0; j < 100; j++)
		{
			for (int k = 0; k < 100; k++)
			{
				dic[j][k] = -1;
			}
		}
	}
	float value(float red, float black)
	{
		float a;
		float b;
		if (red == 0 && black > 0)
		{
			b = black - red;
			return b;
		}
		if (red > 0 && black == 0)
		{
			return 0;
		}
		if (red > 0 && black > 0)
		{
			a = (red / (red + black) * value(red - 1, black)) + (black / (red + black) * value(red, black - 1));
			b = black - red;
			if (a > b)
				return a;
			else
				return b;
		}
	}
	float m_value(float red, float black)
	{
		float a;
		float b;
		if (red == 0 && black > 0)
		{
			b = black - red;
			return b;
		}
		if (red > 0 && black == 0)
		{
			return 0;
		}
		if (red > 0 && black > 0)
		{
			float c, d;
			if (seek(red - 1, black) != -1)
				c = seek(red - 1, black);
			else
			{
				in(red - 1, black,m_value(red - 1, black));
				c = seek(red - 1, black);
			}
			if (seek(red, black - 1) != -1)
				d = seek(red, black - 1);
			else
			{
				in(red, black - 1,m_value(red, black - 1));
				d = seek(red, black - 1);
			}
			a = (red / (red + black) * c) +(black/ (red + black) * d);
			b = black - red;
			if (a > b)
				return a;
			else
				return b;
		}
	}
public:
	// Solves the n-Queens problem by (recursive) backtracking
	void result(float n)
	{
		initialize();
		float red = n / 2;
		float black = n / 2;
		num=m_value(red,black);
		cout << "Value of the game = " << num<< endl;
	}
};
int main(int argc, char** argv)
{
	clock_t start, end;
	Board x;
	int card;
	cout << "Total Number of Cards = ";
	cin >> card;
	start = clock();
	x.result(card);
	end = clock();   //end
	std::cout << "Running time = " << double(end - start) / CLOCKS_PER_SEC << "s" << std::endl;  //print the time
	return 0;
}

