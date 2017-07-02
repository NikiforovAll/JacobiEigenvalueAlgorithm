#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include "JacobiEigenvalueAlgorithm.h"
#include <iostream>
#include "iostream"
#include <string>
#include "string"
#include <boost/numeric/ublas/io.hpp>
#include <iomanip>
#include <omp.h>

void rotateJacobi3(
	boost::numeric::ublas::matrix<double> &S,
	int row, int col)
{
	int q = row;
	int p = col;
	double t, c, s, theta;
	int n = S.size1();
	if (row == col) return;
	double tau = (S(q, q) - S(p, p)) / (2 * S(p, q));
	t = (tau >= 0 ? 1.0 : -1.0) / (abs(tau) + sqrt(1 + tau*tau));
	c = 1.0 / sqrt(1 + t*t);
	s = t*c;
	//0
	boost::numeric::ublas::matrix<double> R(n, n);
	R = boost::numeric::ublas::identity_matrix<double>(n, n);
	R(row, row) = R(col, col) = c; R(row, col) = s; R(col, row) = -s;
	S = prod(S, trans(R));
	S = prod(R, S);
	// 1 
	/*boost::numeric::ublas::matrix<double> A(S.size1(), S.size2());
	for (size_t i = 0; i<A.size1(); i++) {
		for (size_t j = 0; j<A.size2(); j++) {
			A(i, j) = S(i, j);
		}
	}
	for (int i = 0; i < S.size2(); i++) {
		if (i == q || i == p) continue;
		S(i, q) = c*A(i, q) - s*A(i, p);
		S(q, i) = S(i, q);
 	}

	for (int i = 0; i < S.size2(); i++) {
		if (i == q || i == p) continue;
		S(i, p) = c*A(i, p) + s*A(i, q);
		S(i, p) = S(p, i);
	}

	S(q, q) = pow(c, 2)*A(q, q) - 2 * c*s*A(q, p) + pow(s, 2)*A(p, p);
	S(p, p) = pow(s, 2)*A(q, q) + 2 * c*s*A(q, p) + pow(c, 2)*A(p, p);
	S(q, p) = (pow(c, 2) - pow(s, 2))*A(q, p)+s*c*(A(q,q)- A(p,p));
	S(p, q) = S(q, p);*/
	// 2
	//double x, y;
	//// left multiplication
	//for (int i = 0; i < S.size1(); i++) {
	//	if ((p + i - 1) >= S.size1()) break;
	//	x = S(p, p + i - 1);
	//	y = S(q, p + i - 1);
	//	S(p, p + i - 1) = c*x - s*y;
	//	S(q, p + i - 1) = s*x + c*y;
	//}
	//// right multiplication 
	//for (int i = 0; i < S.size2(); i++) {
	//	if ((p + i - 1) >= S.size1()) break;
	//	x = S(p + i - 1, p);
	//	y = S(p + i - 1, q);
	//	S(p + i - 1, p) = c*x - s*y;
	//	S(p + i - 1, q) = s*x + c*y;
	//}

}

int jacm3sync(
	boost::numeric::ublas::matrix<double> &S,
	boost::numeric::ublas::vector<double> &e,
	boost::numeric::ublas::matrix<double>  &U,
	int &iter
) {
	int n = S.size1();
	boost::numeric::ublas::matrix<double> copied(n, n);

	for (size_t i = 0; i<copied.size1(); i++) {
		for (size_t j = 0; j<copied.size2(); j++) {
			copied(i, j) = S(i,j);
		}
	}
	std::cout << "Start matrix:";
	debugMatrixForm(S);

	int p = 3;
	int q = p+1;
		rotateRowCol(S, U, q, p);
		debugMatrixForm(S);
	/*	rotateJacobi3(copied, p, q);
		debugMatrixForm(copied);
	*/
	return 0;
}


