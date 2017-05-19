#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include "JacobiEigenvalueAlgorithm.h"

//LOG:
//   date: [03.05.2017], status: [does not work] if isOptimized = true, but rotation works for JacobiGivensSync.cpp
//   date: [03.05.2017], status : [works] if isOptimized = true, but slow for size(matrix) >= 64

int jacobiSync(
	boost::numeric::ublas::matrix<double> &S,
	boost::numeric::ublas::vector<double> &e,
	boost::numeric::ublas::matrix<double>  &U,
	int &iter,
	bool isOptimized)
{
	iter = 0;
	int col, row;
	bool iterating = true;
	int n = S.size1();
	if (S.size2() != n)
	{
		return -1;
	}
	boost::numeric::ublas::matrix<double> M(n, n);
	e = boost::numeric::ublas::zero_vector<double>(n);
	U = boost::numeric::ublas::identity_matrix<double>(n, n);
	while (iterating)
	{
		M = S;
		abs(M);
		iter++;
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
		if (isOptimized) {
			rotateColRowJacobi(S, row, col);
		}
		else {
			rotateRowCol(S, U, row, col);
		}
		//cout<<row<<" row|col "<<col <<" sum: "<< sumOffDiagonal(S) << endl;
		//if (Smax < _EPS * norm_frobenius(S)) iterating = false;
		if (sumOffDiagonal(S) < _EPS) iterating = false;
	}

	for (int i = 0; i < n; i++) e(i) = S(i, i);

	return 0;
}