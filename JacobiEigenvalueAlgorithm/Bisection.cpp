#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/progress.hpp>
#include "boost/format.hpp"
#include <omp.h>
#include <chrono>
#include <fstream>
#include <sstream>
#include <string>
#include <list>
#include "JacobiEigenvalueAlgorithm.h"

#include <cmath> 

#ifdef omptest
#include <omp.h>
#endif

#include < stdio.h>
#include "f2c.h"
#include "clapack.h"

#include <float.h>                 // required for DBL_EPSILON

#define omptest false

void classic_bisection(
	boost::numeric::ublas::vector<real> diagonal,
	boost::numeric::ublas::vector<real> offdiagonal,
	boost::numeric::ublas::vector<real> &result);
void modified_bisection(
	boost::numeric::ublas::vector<real> diagonal,
	boost::numeric::ublas::vector<real> offdiagonal,
	boost::numeric::ublas::vector<real> &eigenvalues);
//UTIL
void estimate_eigen_interval(
	boost::numeric::ublas::vector<real> diagonal,
	boost::numeric::ublas::vector<real> offdiagonal,
	real& left_boundary,
	real& right_boundary
);
int sturm_sequence(real d[], real off[], real x, int n);

//(bisection)
//http://physics.oregonstate.edu/~landaur/nacphy/lapack/routines/sstebz.html
void bisection_test(
	boost::numeric::ublas::matrix<double> M,
	std::string isWriteToConsole,
	std::ofstream fp_outs[1],
	int i) {
	// PREP
	int iter = 0;
	auto begin = std::chrono::high_resolution_clock::now();
	auto end = std::chrono::high_resolution_clock::now();
	integer matrixSize = M.size1();
	//real *diagonal = new real[matrixSize];
	//real* offdiagonal = new real[matrixSize];
	boost::numeric::ublas::vector<real> result(matrixSize);
	boost::numeric::ublas::vector<real> diagonal(matrixSize);
	boost::numeric::ublas::vector<real> offdiagonal(matrixSize);
	for (int j = 0; j < matrixSize; j++) {
		diagonal[j] = M(j, j);
	}
	offdiagonal[0] = 0.0;
	for (int j = 1; j < matrixSize; j++) {
		offdiagonal[j] = M(j - 1, j);
	}
	//real *result = new real[matrixSize];

	// BEGIN TEST
	begin = std::chrono::high_resolution_clock::now();
	//classic_bisection(diagonal, offdiagonal, result);
	modified_bisection(diagonal, offdiagonal, result);
	end = std::chrono::high_resolution_clock::now();
	// END TEST
	double duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1000000.0;
	// INFO
	if (isWriteToConsole == "true") {
		std::string eig = "[";
		eig += std::to_string(matrixSize);
		eig += "](";
		for (int i = 0; i < matrixSize - 1; i++)
		{
			eig += std::to_string(result[i]);
			eig += ",";
		}
		eig += std::to_string(result[matrixSize - 1]);
		eig += ")";

		writeToAllStreams((boost::format("#%1%: \n") % i).str(), fp_outs);

		writeToAllStreams((boost::format("Name: %1% \nEigenvalues: %2% \nElapsed(ms): %3% \nIter: %4%")
			% "bisection"% eig%duration%iter).str(), fp_outs);

		writeToAllStreams("============================", fp_outs);
	}
}


void modified_bisection(
	boost::numeric::ublas::vector<real> diagonal,
	boost::numeric::ublas::vector<real> offdiagonal,
	boost::numeric::ublas::vector<real> &eigenvalues)
{
	real left_boundary;
	real right_boundary;
	estimate_eigen_interval(diagonal, offdiagonal, left_boundary, right_boundary);

}
////////////////////////////////////////////////////////////////////////////////
// classic_bisection                                                    //
// "Calculation of the Eigenvalues of a Symmetric Tridiagonal Matrix" by the Method of Bisection" W.Barth, R.S. Martin, J.H. Wilkinson. 
////////////////////////////////////////////////////////////////////////////////
//      
void classic_bisection(
	boost::numeric::ublas::vector<real> diagonal,
	boost::numeric::ublas::vector<real> offdiagonal,
	boost::numeric::ublas::vector<real> &eigenvalues)
{
	int n = diagonal.size();
	int i, j, k;
	double tolerance;
	double left_boundary;
	double right_boundary;
	double x;
	double q;
	double xlower, xupper, xmid;
	double *lowerbounds;
	double relative_tolerance = 10e-05;
	double epsilon;
	lowerbounds = (double*)malloc(n * sizeof(double));
	// Use Gerschgorin's Theorem to Find Upper and Lower Bounds for
	// All Eigenvalues.

	offdiagonal[0] = 0.0;
	right_boundary = diagonal[n - 1] + fabs(offdiagonal[n - 1]);
	left_boundary = diagonal[n - 1] - fabs(offdiagonal[n - 1]);
	for (i = n - 2; i >= 0; i--) {
		x = fabs(offdiagonal[i]) + fabs(offdiagonal[i + 1]);
		if ((diagonal[i] + x) > right_boundary) right_boundary = diagonal[i] + x;
		if ((diagonal[i] - x) < left_boundary) left_boundary = diagonal[i] - x;
	}

	//estimate_eigen_interval(diagonal, offdiagonal, left_boundary, right_boundary);
	epsilon = ((right_boundary + left_boundary) > 0.0) ? left_boundary : right_boundary;
	epsilon *= DBL_EPSILON;
	if (relative_tolerance < epsilon) relative_tolerance = epsilon;
	epsilon = 0.5 * relative_tolerance + 7.0 * epsilon;

	for (i = 0; i < n; i++) {
		eigenvalues[i] = right_boundary;
		lowerbounds[i] = left_boundary;
	}
	// Find all eigenvalues from largest to smallest storing
	// from smallest to largest. 

	xupper = right_boundary;
	#ifdef omptest
	#pragma omp parallel for
	#endif
	for (k = n - 1; k >= 0; k--) {
		xlower = left_boundary;
		for (i = k; i >= 0; i--)
			if (xlower < lowerbounds[i]) { xlower = lowerbounds[i]; break; }

		if (xupper > eigenvalues[k]) xupper = eigenvalues[k];

		tolerance = 2.0 * DBL_EPSILON * (fabs(xlower) + fabs(xupper)) + relative_tolerance;

		while ((xupper - xlower) > tolerance) {
			xmid = 0.5 * (xupper + xlower);
			j = sturm_sequence(&diagonal[0], &offdiagonal[0], xmid, n) - 1;
			if (j < k) {
				if (j < 0) { xlower = lowerbounds[0] = xmid; }
				else {
					xlower = lowerbounds[j + 1] = xmid;
					if (eigenvalues[j] > xmid) eigenvalues[j] = xmid;
				}
			}
			else xupper = xmid;
			tolerance = 2.0 * DBL_EPSILON * (fabs(xlower) + fabs(xupper)) + relative_tolerance;
		};
		eigenvalues[k] = 0.5 * (xupper + xlower);
		
	}
}
////////////////////////////////////////////////////////////////////////////////
// estimate_eigen_interval                                                    //
////////////////////////////////////////////////////////////////////////////////
//          
void estimate_eigen_interval(
	//TBD: modify with bnorm max(gl, gu)
	boost::numeric::ublas::vector<real> diagonal,
	boost::numeric::ublas::vector<real> offdiagonal,
	real& left_boundary,
	real& right_boundary
) 
{
	int matrixSize = diagonal.size();
	boost::numeric::ublas::vector<real> radius_i(matrixSize);
	boost::numeric::ublas::vector<real> left_estimation(matrixSize);
	boost::numeric::ublas::vector<real> right_estimation(matrixSize);

	for (size_t i = 0; i < matrixSize; i++)
	{
		if (i != matrixSize - 1)
		{
			radius_i[i] = fabs(offdiagonal[i]) + fabs(offdiagonal[i+1]);
		}
		else
		{
			radius_i[i] = fabs(offdiagonal[i]);
		}
	}

	for (size_t i = 0; i < matrixSize; i++)
	{
		left_estimation[i] = diagonal[i] - radius_i[i];
		right_estimation[i] = diagonal[i] + radius_i[i];
	}


	left_boundary = *std::min_element(std::begin(left_estimation), std::end(left_estimation));
	right_boundary = *std::max_element(std::begin(right_estimation), std::end(right_estimation));
}


////////////////////////////////////////////////////////////////////////////////
//  int Sturm_Sequence                                                        //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
static int sturm_sequence(real d[], real off[], real x, int n)
{
	double q = 1.0;
	int k = 0;
	int i;

	for (i = 0; i < n; i++) {
		if (q == 0.0)
			q = d[i] - x - fabs(off[i]) / DBL_EPSILON;
		else
			q = d[i] - x - off[i] * off[i] / q;
		if (q < 0.0) k++;
	}
	return k;
}

void compute_group_bisect(
	int p,
	int q,
	double left_boundary,
	double right_boundary,
	boost::numeric::ublas::vector<real> &result
)
{

}