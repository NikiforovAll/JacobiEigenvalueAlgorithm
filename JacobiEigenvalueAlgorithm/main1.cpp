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

#include < stdio.h>
#include "f2c.h"
#include "clapack.h"


//using namespace std;
//using namespace boost::numeric::ublas;
//using namespace std::chrono;

// sync Jacobi algorithm: find max O(n2) + full matrix multiplication O(n2)
void JacobiSyncTest(boost::numeric::ublas::matrix<double> M, std::string isWriteToConsole, std::ofstream fp_outs[1], int i) {
	int iter;
	auto begin = std::chrono::high_resolution_clock::now();
	auto end = std::chrono::high_resolution_clock::now();
	boost::numeric::ublas::matrix<double> U(M.size1(), M.size2());
	boost::numeric::ublas::vector<double> e(M.size1());
	double duration = 0;
	//NAME:  jacobiPseudoAsync
	// BEGIN TEST
	begin = std::chrono::high_resolution_clock::now();
	jacobiSync(M, e, U, iter, false);
	end = std::chrono::high_resolution_clock::now();
	// END TEST
	std::sort(e.begin(), e.end(), std::less<double>());
	duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1000000.0;
	if (isWriteToConsole == "true") {
		writeToAllStreams((boost::format("#%1%: \n") % i).str(), fp_outs);

		writeToAllStreams((boost::format("Name: %1% \nEigenvalues: %2% \nElapsed(ms): %3% \nIter: %4%")
			% "JacobiSync"% e%duration%iter).str(), fp_outs);

		writeToAllStreams("============================", fp_outs);
	}
}

// Jacobi algorithm + Givens rotations (matrix clean up)
void JacobiGivensSyncTest(boost::numeric::ublas::matrix<double> M, std::string isWriteToConsole, std::ofstream fp_outs[1], int i) {
	int iter;
	auto begin = std::chrono::high_resolution_clock::now();
	auto end = std::chrono::high_resolution_clock::now();
	boost::numeric::ublas::matrix<double> U(M.size1(), M.size2());
	boost::numeric::ublas::vector<double> e(M.size1());
	double duration = 0;
	// BEGIN TEST
	begin = std::chrono::high_resolution_clock::now();
	jacobiPseudoAsync(M, e, U, iter);
	end = std::chrono::high_resolution_clock::now();
	// END TEST

	duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1000000.0;
	std::sort(e.begin(), e.end(), std::less<double>());
	// INFO
	if (isWriteToConsole == "true") {
		writeToAllStreams((boost::format("#%1%: \n") % i).str(), fp_outs);

		writeToAllStreams((boost::format("Name: %1% \nEigenvalues: %2% \nElapsed(ms): %3% \nIter: %4%")
			% "JacobiGivensSync"% e%duration%iter).str(), fp_outs);

		writeToAllStreams("============================", fp_outs);
	}
}

void JacobiAsyncTest(boost::numeric::ublas::matrix<double> M, std::string isWriteToConsole, std::ofstream fp_outs[1], int i) {
	int iter;
	auto begin = std::chrono::high_resolution_clock::now();
	auto end = std::chrono::high_resolution_clock::now();
	boost::numeric::ublas::matrix<double> U(M.size1(), M.size2());
	boost::numeric::ublas::vector<double> e(M.size1());
	double duration = 0;
	// BEGIN TEST
	begin = std::chrono::high_resolution_clock::now();
	jacobiAsync(M, e, U, iter);
	end = std::chrono::high_resolution_clock::now();
	// END TEST
	duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1000000.0;
	std::sort(e.begin(), e.end(), std::less<double>());
	// INFO
	if (isWriteToConsole == "true") {
		writeToAllStreams((boost::format("#%1%: \n") % i).str(), fp_outs);

		writeToAllStreams((boost::format("Name: %1% \nEigenvalues: %2% \nElapsed(ms): %3% \nIter: %4%")
			% "jabm1async"% e%duration%iter).str(), fp_outs);

		writeToAllStreams("============================", fp_outs);
	}

}

void jacm3syncTest(boost::numeric::ublas::matrix<double> M, std::string isWriteToConsole, std::ofstream fp_outs[1], int i) {
	int iter;
	auto begin = std::chrono::high_resolution_clock::now();
	auto end = std::chrono::high_resolution_clock::now();
	boost::numeric::ublas::matrix<double> U(M.size1(), M.size2());
	boost::numeric::ublas::vector<double> e(M.size1());
	double duration = 0;
	// BEGIN TEST
	begin = std::chrono::high_resolution_clock::now();
	jacm2sync(M, e, U, iter);
	end = std::chrono::high_resolution_clock::now();
	// END TEST

	duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1000000.0;
	std::sort(e.begin(), e.end(), std::less<double>());
	// INFO
	if (isWriteToConsole == "true") {
		writeToAllStreams((boost::format("#%1%: \n") % i).str(), fp_outs);

		writeToAllStreams((boost::format("Name: %1% \nEigenvalues: %2% \nElapsed(ms): %3% \nIter: %4%")
			% "jacm2sync"% e%duration%iter).str(), fp_outs);

		writeToAllStreams("============================", fp_outs);
	}
}

//SSTEQR (QR/QL based)
//http://physics.oregonstate.edu/~landaur/nacphy/lapack/routines/ssteqr.html
void ssteqr_lapacktest(boost::numeric::ublas::matrix<double> M, std::string isWriteToConsole, std::ofstream fp_outs[1], int i) {
	// PREP
	int iter = 0;
	auto begin = std::chrono::high_resolution_clock::now();
	auto end = std::chrono::high_resolution_clock::now();
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
	begin = std::chrono::high_resolution_clock::now();
	ssteqr_(&c4, &currMatrixSize, diagonal, offdiagonal, dummy, &currMatrixSize, dummy, &info);
	end = std::chrono::high_resolution_clock::now();
	// END TEST
	double duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1000000.0;
	// INFO
	if (isWriteToConsole == "true") {
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

		writeToAllStreams((boost::format("#%1%: \n") % i).str(), fp_outs);

		writeToAllStreams((boost::format("Name: %1% \nEigenvalues: %2% \nElapsed(ms): %3% \nIter: %4%")
			% "ssteqr" % eig % duration % iter).str(), fp_outs);

		writeToAllStreams("============================", fp_outs);
	}
	delete[] diagonal;
	delete[] offdiagonal;
	delete[] dummy;
}
//SSTEBZ (bisection)
//http://physics.oregonstate.edu/~landaur/nacphy/lapack/routines/sstebz.html
void sstebz_lapacktest(boost::numeric::ublas::matrix<double> M, std::string isWriteToConsole, std::ofstream fp_outs[1], int i) {
	// PREP
	int iter = 0;
	auto begin = std::chrono::high_resolution_clock::now();
	auto end = std::chrono::high_resolution_clock::now();
	integer matrixSize = M.size1();
	real *diagonal = new real[matrixSize];
	real* offdiagonal = new real[matrixSize];
	for (int j = 0; j < matrixSize; j++) {
		diagonal[j] = M(j, j);
	}
	for (int j = 0; j < matrixSize - 1; j++) {
		offdiagonal[j] = M(j + 1, j);
	}
	integer info = 0;
	real abstol = 0;
	char c4, c5;
	c4 = 'A';
	c5 = 'E';
	real dummy1 = 0;
	integer dummy2 = 0;
	integer numberOfEigenvaluesFound = 0, numberOfdiagonalBlocks = 0;
	real *result = new real[matrixSize];
	integer *iblock = new integer[matrixSize];
	integer *split = new integer[matrixSize];
	real *workspace = new real[matrixSize * 4];
	integer *iworkspace = new integer[matrixSize * 3];

	// BEGIN TEST
	begin = std::chrono::high_resolution_clock::now();
	sstebz_(&c4, &c5, &matrixSize,
		&dummy1, &dummy1, &dummy2, &dummy2,
		&abstol, diagonal, offdiagonal,
		&numberOfEigenvaluesFound, &numberOfdiagonalBlocks,
		result, iblock, split, workspace, iworkspace, &info);
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
			% "sstebz"% eig%duration%iter).str(), fp_outs);

		writeToAllStreams("============================", fp_outs);
	}
	delete[] diagonal, offdiagonal, dummy1, result, iblock, split, workspace, iworkspace;

}

int main2(int argc, char **argv)
{

	//TBD: review code for memory overhead
	int startIndex = 1;
	int shift = 1;

	int numberOfMatrix = startIndex + shift;
	std::string isWriteToConsole = "true";
	if (argc > 1 && argv) {

		numberOfMatrix = std::stoi(argv[1]);
		isWriteToConsole = argv[2];
	}

	std::ofstream fp_outs[1];
	fp_outs[0].open("output.txt", std::ios::out);
	boost::numeric::ublas::matrix<double>*MatrixArray = readFromSample(numberOfMatrix, "input.txt");
	std::cout << "INFO: read completed." << std::endl;

	for (int i = startIndex; i < startIndex + shift; i++)
	{
		//JacobiSyncTest(MatrixArray[i], isWriteToConsole, fp_outs, i);
		//JacobiGivensSyncTest(MatrixArray[i], isWriteToConsole, fp_outs, i);
		//ssteqr_lapacktest(MatrixArray[i], isWriteToConsole, fp_outs, i);
		jacm3syncTest(MatrixArray[i], isWriteToConsole, fp_outs, i);
		//sstebz_lapacktest(MatrixArray[i], isWriteToConsole, fp_outs, i);
		//doesn't work 
		//JacobiAsyncTest(MatrixArray[i], isWriteToConsole, fp_outs, i);
	}
	std::cout << "INFO: work finished." << std::endl;
	std::cin.get();
	//return exit_success;
	return 0;

}
