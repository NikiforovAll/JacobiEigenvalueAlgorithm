#include <boost/numeric/ublas/matrix.hpp>
#include <string>
//using namespace std;
//using namespace boost::numeric;
//using namespace boost;

//dirty approach 

//Util 
boost::numeric::ublas::matrix<double>* readFromSample(int num, std::string filename);
void writeToAllStreams(std::string str, std::ofstream stream[1]);
void abs(boost::numeric::ublas::matrix<double> &M);
double sumOffDiagonal(boost::numeric::ublas::matrix<double> &S);
void findMax(boost::numeric::ublas::matrix<double> &M, int &row, int &col);

//Methods

//JacobiGivensSync.cpp
int jacobiPseudoAsync(
	boost::numeric::ublas::matrix<double> &S,
	boost::numeric::ublas::vector<double> &e, 
	boost::numeric::ublas::matrix<double>  &U,
	int &iter
);

int jacobiSync(
	boost::numeric::ublas::matrix<double> &S,
	boost::numeric::ublas::vector<double> &e,
	boost::numeric::ublas::matrix<double>  &U,
	int &iter,
	bool isOptimized
);

int jacobiAsync(
	boost::numeric::ublas::matrix<double> &S, 
	boost::numeric::ublas::vector<double> &e, 
	boost::numeric::ublas::matrix<double>  &U, int &iter);

//Rotations 
void rotateColRowJacobi(
	boost::numeric::ublas::matrix<double> &S, 
	int row, 
	int col);
void rotateColRowGivens(
	boost::numeric::ublas::matrix<double> &S, 
	int row, 
	int col, 
	int len);
void rotateRowCol(
	boost::numeric::ublas::matrix<double> &S, 
	boost::numeric::ublas::matrix<double> &U, int row, int col);
#define _EPS 10e-03
