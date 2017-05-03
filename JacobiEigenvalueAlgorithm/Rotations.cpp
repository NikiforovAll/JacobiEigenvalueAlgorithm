
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
using namespace std;
using namespace boost::numeric::ublas;
using namespace boost;
#include "JacobiEigenvalueAlgorithm.h"

// simple rotation 
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

void rotateColRowJacobi(matrix<double> &S, int row, int col) {
	int p = col;
	int q = row;
	double tau = (S(q, q) - S(p, p)) / (2 * S(p, q));
	double c, s, t;
	if (tau == 0) {
		c = 1 / sqrt(2);
		s = c;
	}
	else {
		t = (tau >= 0 ? 1.0 : -1.0) / (abs(tau) + sqrt(1 + tau*tau));
		c = 1.0 / sqrt(1 + t*t);
		s = t*c;
	}
	double x, y;
	// left multiplication
	for (int i = 1; i <= 3; i++) {
		if ((p + i - 1) >= S.size1()) break;
		x = S(p, p + i - 1);
		y = S(q, p + i - 1);
		S(p, p + i - 1) = c*x - s*y;
		S(q, p + i - 1) = s*x + c*y;
	}
	// right multiplication 
	for (int i = 1; i <= 3; i++) {
		if ((p + i - 1) >= S.size1()) break;
		x = S(p + i - 1, p);
		y = S(p + i - 1, q);
		S(p + i - 1, p) = c*x - s*y;
		S(p + i - 1, q) = s*x + c*y;
	}
}

void rotateColRowGivens(matrix<double> &S, int row, int col, int len) {
	double x, y, t, s, c;
	int p, q;
	p = row;
	q = col;
	x = S(p, p - 1);
	y = S(q, p - 1);
	t = sqrt(x*x + y*y);
	s = -y / t;
	c = x / t;
	S(p, p - 1) = t;
	S(q, p - 1) = 0;
	// left multiplication 
	for (int i = 1; i <= len; i++) {
		x = S(p, p + i - 1);
		y = S(q, p + i - 1);
		S(p, p + i - 1) = c*x - s*y;
		S(q, p + i - 1) = s*x + c*y;
	}
	S(p - 1, p) = t;
	S(p - 1, q) = 0;
	// right multiplication 
	for (int i = 1; i <= len; i++) {
		x = S(p + i - 1, p);
		y = S(p + i - 1, q);
		S(p + i - 1, p) = c*x - s*y;
		S(p + i - 1, q) = s*x + c*y;
	}
}