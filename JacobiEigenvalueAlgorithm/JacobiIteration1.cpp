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
#include < stdio.h>
#include "f2c.h"
#include "clapack.h"

#define _EPS 10e-03
void ssteqr_debug(boost::numeric::ublas::matrix<double> M) {
	// PREP
	integer currMatrixSize = M.size1();
	real *diagonal = new real[currMatrixSize];
	real* offdiagonal = new real[currMatrixSize];
	for (int j = 0; j < currMatrixSize; j++) {
		diagonal[j] = M(j, j);
	}
	for (int j = 0; j < currMatrixSize - 1; j++) {
		offdiagonal[j] = M(j + 1, j);
	}
	real* dummy = new real[currMatrixSize * 2];
	integer info = 0;
	char c4;
	c4 = 'N';
	// BEGIN TEST
	// END TEST
	// INFO
	std::string eig = "[";
	eig += std::to_string(currMatrixSize);
	eig += "](";
	for (int i = 0; i < currMatrixSize - 1; i++)
	{
		eig += std::to_string(diagonal[i]);
		eig += ",";
	}
	eig += std::to_string(diagonal[currMatrixSize - 1]);
	eig += ")";

	std::cout << eig;
	delete[] diagonal;
	delete[] offdiagonal;
	delete[] dummy;
}
int rotateJacobi2(boost::numeric::ublas::matrix<double> &S, int row, int col) {
	int q = row;
	int p = col;
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
	int n = S.size1();
	boost::numeric::ublas::matrix<double> R(S.size1(), S.size1());
	R = boost::numeric::ublas::identity_matrix<double>(n, n);
	R(row, row) = R(col, col) = c; R(row, col) = s; R(col, row) = -s;
	S = prod(S, trans(R));
	S = prod(R, S);
	return 0;
}


int jacm2sync(
	boost::numeric::ublas::matrix<double> &S,
	boost::numeric::ublas::vector<double> &e,
	boost::numeric::ublas::matrix<double>  &U,
	int &iter
) {
	int n = S.size1();
	iter = 0;
	//omp_set_num_threads(4);

	int thread_count = omp_get_num_threads();
	bool iterating = true;

	while (iterating)
	{
		debugMatrixForm(S);
		ssteqr_debug(S);

		//#pragma omp parallel for shared(S) 
		for (int pp = 0; pp < n - 2; pp += 2) {
			int p = pp;
			int q = p + 1;
			////printf("In thread %d p = %d q = %d\n", omp_get_thread_num(), p, q);
			//std::cout << "S(p, q)" << S(p, q) << " | " << "S(q, p)" << S(q, p) << std::endl;
			//debugMatrixForm(S);
			//std::cout << p <<" p|q "<< q << std::endl;

			rotateJacobi2(S, q, p);
			debugMatrixForm(S);
			ssteqr_debug(S);

			//std::cout << "S(p, q)" << S(p, q) << " | "<< "S(q, p)" << S(q, p) << std::endl;
			//debugMatrixForm(S);
			//std::cout << p << " p|q " << q << std::endl;
			//std::cout << "====" << std::endl;
		}

		//#pragma omp parallel for shared(S)
		for (int pp = 1; pp < n - 1; pp += 2) {
			int p = pp;
			int q = p + 1;
			/*std::cout << "S(p, q)" << S(p, q) << " | " << "S(q, p)" << S(q, p) << std::endl;
			debugMatrixForm(S);
			std::cout << p <<" p|q "<< q << std::endl;*/
			
			
			rotateJacobi2(S, q, p);

		/*	std::cout << "S(p, q)" << S(p, q) << " | "<< "S(q, p)" << S(q, p) << std::endl;
			debugMatrixForm(S);
			std::cout << p << " p|q " << q << std::endl;
			std::cout << "====" << std::endl;*/
		}

		debugMatrixForm(S);

		int q = 0;
		int p = 0;
		for (int p = 1; p < n - 2; p++) {
			q = p + 1;
			if (abs(S(q, p - 1)) > _EPS*sqrt(S(q, q)*S(q, q) + S(p, p)*(p, p))) {
				rotateColRowGivens(S, p, q, 3);
				//debugMatrixForm(S);
			}
		}
		p = n - 2;
		q = p + 1;
		//last step
		if (abs(S(p, q))> _EPS* sqrt(S(q, q)*S(q, q) + S(p, p)*(p, p))) {
			rotateColRowGivens(S, p, q, 2);
			//debugMatrixForm(S);
		}

		debugMatrixForm(S);
		double sum = summOffDiagonal2(S);
		iter++;
		if (sum < _EPS) iterating = false;
	}
	for (int i = 0; i < S.size1(); i++)
	{
		e(i) = S(i, i);
	}
	return 0;
}


