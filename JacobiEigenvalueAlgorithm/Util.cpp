#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/progress.hpp>
#include "boost/format.hpp"
#include <fstream>
using namespace std;
using namespace boost::numeric::ublas;
using namespace boost;

#define _EPS 10e-03
//format: [size, size]((n1,...,size),(n1,...size)...(n1,...size))
// number format: ###.###
matrix<double>* readFromSample(int num, string filename) {
	string line;
	ifstream myfile(filename);
	int count = 0;
	string matrixStr;
	matrix<double>* result;
	result = new (nothrow) matrix<double>[num];
	bool isMatrix = false;
	if (myfile.is_open())
	{
		while (getline(myfile, line))
		{
			if (count < num)
				if (line == "#") {
					bool sizeState = true;
					bool isReadState = false;
					string sizeStr;
					string rowStr = "";
					int size[2];
					int currIndex[2] = { 0, 0 };
					for (int i = 1; i < matrixStr.size() - 1; i++) {
						if (sizeState) {
							if (matrixStr[i] == ']') {
								std::istringstream ss(sizeStr);
								std::string token;
								int currSize = 0;
								while (std::getline(ss, token, ',')) {
									size[currSize] = stoi(token);
									currSize++;
								}
								sizeState = false;
								result[count] = matrix<double>(size[0], size[1]);
								i++;
							}
							else {
								sizeStr += matrixStr[i];
							}
						}
						else {
							if (matrixStr[i] == '(') {
								isReadState = true;
								currIndex[1] = 0;
								continue;
							}
							if (isReadState && matrixStr[i] != ')') {
								rowStr += matrixStr[i];
							}
							else if (matrixStr[i] == ')') {
								std::istringstream ss(rowStr);
								std::string token;
								int currSize = 0;
								while (std::getline(ss, token, ',')) {
									result[count](currIndex[0], currSize) = stod(token);
									currSize++;
								}
								currIndex[0]++;
								isReadState = false;
								rowStr = "";
								i++;
							}
						}
					}
					matrixStr = "";
					count++;
					if (count > num) {
						break;
					}
				}
				else {
					matrixStr += line;
				}

		}
		myfile.close();

	}
	else cout << "Unable to open file";
	return result;
}

//writes to console and specified stream
void writeToAllStreams(string str, ofstream stream[1]) {
	cout << str << endl;
	stream[0] << str << endl;
}

void abs(matrix<double> &M)
{
	int n = M.size1();

	for (int k = 0; k < n; k++)
	{
		for (int i = 0; i < n; i++)
		{
			M(k, i) = abs(M(k, i));
		}
	}
}

// zero diagonal matrix expected 
double sumOffDiagonal(matrix<double> &S) {
	double sum = 0;
	for (int i = 0; i < S.size1(); i++) {
		for (int j = i + 1; j < S.size2() - 1; j++)
		{
			sum += abs(S(i, j));
		}
	}
	return sum;
}

// finding pivot (off diagonal element)
void findMax(matrix<double> &M, int &row, int &col)
{
	double m = M(0, 0);
	row = col = 0;
	int n = M.size1();

	for (int k = 0; k < n; k++)
	{
		for (int l = k; l < n; l++)
		{
			if (M(k, l) > m)
			{
				row = k;
				col = l;
				m = M(k, l);
			}
		}
	}
}