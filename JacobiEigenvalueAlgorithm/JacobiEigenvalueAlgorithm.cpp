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

#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>
#include <sstream>
#include <iterator>
#include <cstdlib>
#include <cmath>

#include "matrix.h"
//#include "gaussian_elimination.h"
#include "parallel_jacobi.h"
#include <vector>

static const float epsilon = 1e-1f;

void find_eigenvalues_parallel_2(matrix* A, std::vector<float>& eigenvalues)
{
	//timer t("run");
	parallel_jacobi::converge_off_threshold sc(1e-6, *A);
	parallel_jacobi::music_permutation pe(A->size());
	parallel_jacobi::run(*A, sc, pe);

	// Eigenvalues are left on the diagonal
	for (int i = 0; i<A->actual_size(); ++i)
		eigenvalues.push_back(A->get(i, i));
}

void make_matrix(matrix& A, float values[]) {
	for (int i = 0; i<A.actual_size(); ++i)
		for (int j = 0; j<A.actual_size(); ++j)
			A(i, j) = values[i*A.actual_size() + j];
}

int main(int argc, char **argv)
{

	float mat[] = { 2.0f, 3.0f, 4.0f, 3.0f, 6.0f, -5.0f, 4.0f, -5.0f, 7.0f };
	matrix A(3);
	make_matrix(A, mat);
	std::vector<float> eigs;
	find_eigenvalues_parallel_2(&A, eigs);
	std::cout << A;
	std::sort(eigs.begin(), eigs.end());
	return 0;
	
}
