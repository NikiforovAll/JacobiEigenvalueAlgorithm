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

#define _EPS 10e-03


int jacm2sync(
	boost::numeric::ublas::matrix<double> &S,
	boost::numeric::ublas::vector<double> &e,
	boost::numeric::ublas::matrix<double>  &U,
	int &iter
) {
	int n = S.size1();
	iter = 0;
	omp_set_num_threads(4);

	int thread_count = omp_get_num_threads();
	bool iterating = true;
	while (iterating)
	{
		#pragma omp parallel for shared(S) 
		for (int pp = 0; pp < n - 2; pp += 2) {
			int p = pp;
			int q = p + 1;
			//printf("In thread %d p = %d q = %d\n", omp_get_thread_num(), p, q);
			/*std::cout << "S(p, q)" << S(p, q) << " | " << "S(q, p)" << S(q, p) << std::endl;
			debugMatrixForm(S);
			std::cout << p <<" p|q "<< q << std::endl;*/

			rotateColRowJacobi(S, q, p);

			/*std::cout << "S(p, q)" << S(p, q) << " | "<< "S(q, p)" << S(q, p) << std::endl;
			debugMatrixForm(S);
			std::cout << p << " p|q " << q << std::endl;
			std::cout << "====" << std::endl;*/
		}

		#pragma omp parallel for shared(S)
		for (int pp = 1; pp < n - 1; pp += 2) {
			int p = pp;
			int q = p + 1;
			/*std::cout << "S(p, q)" << S(p, q) << " | " << "S(q, p)" << S(q, p) << std::endl;
			debugMatrixForm(S);
			std::cout << p <<" p|q "<< q << std::endl;
			*/
			
			rotateColRowJacobi(S, q, p);

			/*std::cout << "S(p, q)" << S(p, q) << " | "<< "S(q, p)" << S(q, p) << std::endl;
			debugMatrixForm(S);
			std::cout << p << " p|q " << q << std::endl;
			std::cout << "====" << std::endl;*/
		}
		//debugMatrixForm(S);
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


	/*	findMax2(S, p, q);
		if (p != q)
		{
			rotateColRowJacobi(S, p, q);
		}*/

		//debugMatrixForm(S);
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

