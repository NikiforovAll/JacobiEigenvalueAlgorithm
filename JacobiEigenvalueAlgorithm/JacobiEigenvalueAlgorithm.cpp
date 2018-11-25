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
#include "parallel_jacobi.h"
#include "timer.h"
#include <vector>

void make_matrix(matrix& A, boost::numeric::ublas::matrix<float> M) {
	for (int i = 0; i<A.actual_size(); ++i)
		for (int j = 0; j<A.actual_size(); ++j)
			A(i, j) = M(i,j);
}

void init_matrix(matrix** mat, int n)
{
	*mat = new matrix(n);
}

// TEST PART

void find_eigenvalues_parallel_jacob_music(matrix* A, std::vector<float>& eigenvalues, int &iter)
{
	double offThreshold = 10e-3;
	parallel_jacobi::converge_off_threshold sc(offThreshold, *A);
	parallel_jacobi::music_permutation pe(A->size());
	parallel_jacobi::run(*A, sc, pe, iter);
	// Eigenvalues are left on the diagonal
	for (int i = 0; i<A->actual_size(); ++i)
		eigenvalues.push_back(A->get(i, i));
}

void parallel_jacob_musictest(boost::numeric::ublas::matrix<float> M, std::string isWriteToConsole, std::ofstream fp_outs[1], int i) {
	int iter;
	auto begin = std::chrono::high_resolution_clock::now();
	auto end = std::chrono::high_resolution_clock::now();

	double duration = 0;
	matrix* A = 0; 
	init_matrix(&A, M.size1());
	make_matrix(*A, M);
	std::vector<float> e;
	
	// BEGIN TEST
	begin = std::chrono::high_resolution_clock::now();
	find_eigenvalues_parallel_jacob_music(A, e, iter);
	end = std::chrono::high_resolution_clock::now();
	// END TEST
	std::sort(e.begin(), e.end());
	duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1000000.0;
	boost::numeric::ublas::vector<float> eigs(M.size1());
	for (size_t i = 0; i < e.size(); i++)
	{
		eigs(i) = e[i];
	}
	
	// INFO
	if (isWriteToConsole == "true") {
		std::string eig = "[";
		eig += std::to_string(M.size1());
		eig += "](";
		for (int i = 0; i < M.size1() - 1; i++)
		{
			eig += std::to_string(e[i]);
			eig += ",";
		}
		eig += std::to_string(e[M.size1() - 1]);
		eig += ")";

		writeToAllStreams((boost::format("#%1%: \n") % i).str(), fp_outs);

		writeToAllStreams((boost::format("Name: %1% \nEigenvalues: %2% \nElapsed(ms): %3% \nIter: %4%")
			% "parallel_jacob_music"% eig%duration%iter).str(), fp_outs);

		writeToAllStreams("============================", fp_outs);
	}
}




//int main(int argc, char **argv)
//{
//	int startIndex = 0;
//	int shift = 3; // max 0 - 6 - 512
//	int numberOfMatrix = startIndex + shift;
//	std::string isWriteToConsole = "true";
//	if (argc > 1 && argv) {
//
//		numberOfMatrix = std::stoi(argv[1]);
//		isWriteToConsole = argv[2];
//	}
//
//	//omp_set_num_threads(8);
//
//	std::ofstream fp_outs[1];
//	std::string prefix = "1";
//	fp_outs[0].open("output"+ prefix + ".txt", std::ios::out);
//	boost::numeric::ublas::matrix<double>*MatrixArray = readFromSample(numberOfMatrix, "input.txt");
//	std::cout << "INFO: read completed." << std::endl;
//	std::cout << "INFO: ";
//	#ifdef omptest
//		std::cout << "\nRunning parallel jacobi on " << omp_get_max_threads()
//			<< " threads.\n";
//	#else
//		std::cout << "\nRunning serial jacobi.\n";
//	#endif
//	omp_set_num_threads(8);
//	std::cout.precision(12);
//	for (int i = startIndex; i < startIndex + shift; i++)
//	{
//		// QR 
//		ssteqr_lapacktest(MatrixArray[i], isWriteToConsole, fp_outs, i);
//		// Divide&Conquer
//		stedc_lapacktest(MatrixArray[i], isWriteToConsole, fp_outs, i);
//
//		parallel_jacob_musictest(MatrixArray[i], isWriteToConsole, fp_outs, i);
//	}
//	return 0;	
//}

int main(int argc, char **argv)
{
//	int startindex = 2;
//	int shift = 5; // max 0 - 6 - 512
	int startindex = 4;
	int shift = 1;
	int numberofmatrix = startindex + shift;
	std::string IsWriteToConsole = "true";
	if (argc > 1 && argv) {

		numberofmatrix = std::stoi(argv[1]);
		IsWriteToConsole = argv[2];
	}

	int nthreads = 8;
	//omp_set_num_threads(nthreads);
	omp_set_num_threads(nthreads);
	std::ofstream fp_outs[1];
	
	std::string filename = "output"+ std::to_string(nthreads)+".txt";
	//std::string filename = "sstebz_lapacktest.csv";


	fp_outs[0].open(filename, std::ios::out);
	boost::numeric::ublas::matrix<double>*MatrixArray = readFromSample(numberofmatrix, "input.txt");
	//std::cout << "info: read completed." << std::endl;
	//std::cout << "info: ";
	std::cout.precision(12);
	for (int i = startindex; i < startindex + shift; i++)
	{
		//sstebz_lapacktest(MatrixArray[i], IsWriteToConsole, fp_outs, i);
		ssteqr_lapacktest(MatrixArray[i], IsWriteToConsole, fp_outs, i);
		stedc_lapacktest(MatrixArray[i], IsWriteToConsole, fp_outs, i);
		//bisection_test(nthreads, MatrixArray[i], IsWriteToConsole, fp_outs, i);
	}
	return 0;	
}
