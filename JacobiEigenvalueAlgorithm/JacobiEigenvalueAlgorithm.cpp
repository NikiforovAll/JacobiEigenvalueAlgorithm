

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <omp.h>
#include <fstream>
#include <string>
#include <sstream>
#include <list>
#include <boost/progress.hpp>
using namespace std;
using namespace boost::numeric::ublas;
using namespace boost;


#define _EPS 10e-08


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

// finding pivot 
void findMax(matrix<double> &M, int &row, int &col)
{
	double m = M(0, 0);
	row = col = 0;
	int n = M.size1();

	for (int k = 0; k < n; k++)
	{
		for (int l = 0; l < n; l++)
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

// rotation 
void rotateRowCol(matrix<double> &S, matrix<double> &U, int row, int col)
{
	double t, c, s, theta;
	int n = S.size1();
	if (row == col) return;
	theta = (S(row, row) - S(col, col)) / (2 * S(row, col));
	if (theta < _EPS)
	{
		t = 1 / (abs(theta) + sqrt(theta*theta + 1));
	}
	else {
		t = 1 / abs(2 * theta);
	}
	if (theta < 0) t = -t;
	c = 1 / sqrt(t*t + 1);
	s = c*t;
	matrix<double> R(n, n);
	R = identity_matrix<double>(n, n);
	R(row, row) = R(col, col) = c; R(row, col) = s; R(col, row) = -s;
	S = prod(S, trans(R));
	S = prod(R, S);
	U = prod(U, R);
}


int jacobiSync(matrix<double> &S, boost::numeric::ublas::vector<double> &e, matrix<double>  &U, int &iter )
{
	iter = 0;
	int col, row;
	bool iterating = true;
	int n = S.size1();
	if (S.size2() != n)
	{
		return -1;
	}
	matrix<double> M(n, n);
	U = identity_matrix<double>(n, n);
	while (iterating)
	{
		iter++;
		M = S;
		abs(M);
		for (int k = 0; k < n; k++)
		{
			M(k, k) = 0;
		}
		findMax(M, row, col);
		if (row == col)
		{
			for (int i = 0; i < n; i++) e(i) = S(i, i);
			return 0;
		}
		double Smax = S(row, col);
		rotateRowCol(S, U, row, col);
		if (Smax < _EPS * norm_frobenius(S)) iterating = false;
	}

	for (int i = 0; i < n; i++) e(i) = S(i, i);

	return 0;
}



matrix<double>* readFromSample(int num) {
	string line;
	ifstream myfile("sample.txt");
	int count = 0;
	string matrixStr;
	matrix<double>* result;
	result = new (nothrow) matrix<double>[num];
	bool isMatrix = false;
	if (myfile.is_open())
	{
		while (getline(myfile, line))
		{
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

int main(int argc, char **argv)
{
	/*int test(999);

	omp_set_num_threads(2);
#pragma omp parallel reduction(+:test)
	{
#pragma omp critical
		std::cout << "test = " << test << std::endl;
	}
	getchar();*/
	int const numberOfMatrix = 2;

	matrix<double>*MatrixArray = readFromSample(numberOfMatrix);
	for (int i = 0; i < numberOfMatrix; i++)
	{
		boost::timer t;
		int iter; 
		std::cout << "============================" << std::endl;
		matrix<double> M = MatrixArray[i];
		std::cout << "A" << i << std::endl << M << std::endl;
		matrix<double> U(M.size1(), M.size2());
		boost::numeric::ublas::vector<double> e(M.size1());
		jacobiSync(M, e, U, iter);
		std::cout << "Eigenvalues: " << std::endl << e << std::endl;
		std::cout << "U: " << std::endl << U << std::endl;
		std::cout << "Iter: " << std::endl << iter << std::endl;
		std::cout << "Elapsed: "<< std::fixed  << t.elapsed() << std::endl;
	}
	//TBD: output to file ; async solution ; higher precision timer ; 
	return EXIT_SUCCESS;
}