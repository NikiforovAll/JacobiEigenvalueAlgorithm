#include <boost/numeric/ublas/matrix.hpp>
#include <string>
using namespace std;
using namespace boost::numeric::ublas;
using namespace boost;

//dirty approach 

//Util 
matrix<double>* readFromSample(int num, string filename);
void writeToAllStreams(string str, ofstream stream[1]);
void abs(matrix<double> &M);
double sumOffDiagonal(matrix<double> &S);
void findMax(matrix<double> &M, int &row, int &col);

//Methods

//JacobiGivensSync.cpp
int jacobiPseudoAsync(
	matrix<double> &S,
	boost::numeric::ublas::vector<double> &e, matrix<double>  &U,
	int &iter
);

int jacobiSync(
	matrix<double> &S,
	boost::numeric::ublas::vector<double> &e,
	matrix<double>  &U,
	int &iter,
	bool isOptimized
);

int jacobiAsync(matrix<double> &S, boost::numeric::ublas::vector<double> &e, matrix<double>  &U, int &iter);

//Rotations 
void rotateColRowJacobi(matrix<double> &S, int row, int col);
void rotateColRowGivens(matrix<double> &S, int row, int col, int len);
void rotateRowCol(matrix<double> &S, matrix<double> &U, int row, int col);
#define _EPS 10e-03
