#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include "JacobiEigenvalueAlgorithm.h"

#define _EPS 10e-04

//LOG:
//date: [03.05.2017] status: [works]
//date: [03.05.2017] status: [works, but with bad precision]



int jacobiPseudoAsync(
		boost::numeric::ublas::matrix<double> &S,
		boost::numeric::ublas::vector<double> &e, 
		boost::numeric::ublas::matrix<double>  &U,
		int &iter
	) 
{
	int n = S.size1();
	iter = 0;
	//Jacobi rotation
	for (int pp = 0; pp < n - 2; pp++) {
		int p = pp;
		int q = p + 1;
		while (abs(S(p, q)) > _EPS*sqrt(S(q, q)*S(q, q) + S(p, p)*(p, p))) {
			iter++;
			rotateColRowJacobi(S, q, p);
			/*cout << S(p, q)<<endl<<S(q,p)<<endl;
			cout << S << endl;
			cout << p <<" p|q "<< q << endl;*/
			// removing redundant accumulation via Givens rotation
			for (int p = 1; p < n - 2; p++) {
				q = p + 1;
				if (abs(S(q, p - 1)) > _EPS*sqrt(S(q, q)*S(q, q) + S(p, p)*(p, p))) {
					rotateColRowGivens(S, p, q, 3);
				}
			}
			p = n - 2;
			q = p + 1;
			//last step
			if (abs(S(p, q)) > _EPS*sqrt(S(q, q)*S(q, q) + S(p, p)*(p, p))) {
				rotateColRowGivens(S, p, q, 2);
			}
			p = pp;
			q = p + 1;
			//cout <<"iter: "<<endl<<S;
			//cin.get();
		}
	}
	for (int i = 0; i < S.size1(); i++)
	{
		e(i) = S(i, i);
	}
	return 0;
}

