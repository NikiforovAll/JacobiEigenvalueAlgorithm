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

void JacobiGivensSyncTest(boost::numeric::ublas::matrix<double> M, std::string isWriteToConsole, std::ofstream fp_outs[1], int i) {
	int iter;
	auto begin = std::chrono::high_resolution_clock::now();
	auto end = std::chrono::high_resolution_clock::now();
	if (isWriteToConsole == "true") {
		writeToAllStreams((boost::format("#%1%: \n") % i).str(), fp_outs);
		//writeToAllStreams((boost::format("A%1%: \n %2%") % i %M).str(), fp_outs);
	}
	else {
		writeToAllStreams(std::to_string(M.size1()) + "x" + std::to_string(M.size2()), fp_outs);
	}
	boost::numeric::ublas::matrix<double> U(M.size1(), M.size2());
	boost::numeric::ublas::vector<double> e(M.size1());
	double duration = 0;
	//NAME:  jacobiPseudoAsync
	writeToAllStreams("\nJacobiGivensSyncTest", fp_outs);
	// BEGIN TEST
	begin = std::chrono::high_resolution_clock::now();
	jacobiPseudoAsync(M, e, U, iter);
	end = std::chrono::high_resolution_clock::now();
	// END TEST

	duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1000000.0;
	std::sort(e.begin(), e.end(), std::less<double>());
	// INFO
	if (isWriteToConsole == "true") {
		writeToAllStreams((boost::format("Eigenvalues: %1% \nIter %2%\n Elapsed: %3%")
			% e%iter%duration).str(), fp_outs);
	}
	else {
		writeToAllStreams("Elapsed " + std::to_string(duration), fp_outs);
	}

	writeToAllStreams("============================", fp_outs);
}

void JacobiSyncTest(boost::numeric::ublas::matrix<double> M, std::string isWriteToConsole, std::ofstream fp_outs[1], int i) {
	int iter;
	auto begin = std::chrono::high_resolution_clock::now();
	auto end = std::chrono::high_resolution_clock::now();
	if (isWriteToConsole == "true") {
		writeToAllStreams((boost::format("#%1%: \n") % i).str(), fp_outs);
		//writeToAllStreams((boost::format("A%1%: \n %2%") % i %M).str(), fp_outs);
	}
	else {
		writeToAllStreams(std::to_string(M.size1()) + "x" + std::to_string(M.size2()), fp_outs);
	}
	boost::numeric::ublas::matrix<double> U(M.size1(), M.size2());
	boost::numeric::ublas::vector<double> e(M.size1());
	double duration = 0;
	//NAME:  jacobiPseudoAsync
	writeToAllStreams("\nJacobiSyncTest", fp_outs);
	// BEGIN TEST
	begin = std::chrono::high_resolution_clock::now();
	jacobiSync(M, e, U, iter, false);
	end = std::chrono::high_resolution_clock::now();
	// END TEST

	duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1000000.0;

	// INFO
	if (isWriteToConsole == "true") {
		writeToAllStreams((boost::format("Eigenvalues: %1% \n U: %2% \nIter %3%\n Elapsed: %4%")
			% e %U%iter%duration).str(), fp_outs);
	}
	else {
		writeToAllStreams("Elapsed " + std::to_string(duration), fp_outs);
	}

	writeToAllStreams("============================", fp_outs);
}
//
//void JacobiAsyncTest(matrix<double> M, string isWriteToConsole, ofstream fp_outs[1], int i) {
//	int iter;
//	auto begin = std::chrono::high_resolution_clock::now();
//	auto end = std::chrono::high_resolution_clock::now();
//	if (isWriteToConsole == "true") {
//		writeToAllStreams((boost::format("A%1%: \n %2%") % i %M).str(), fp_outs);
//	}
//	else {
//		writeToAllStreams(std::to_string(M.size1()) + "x" + std::to_string(M.size2()), fp_outs);
//	}
//	matrix<double> U(M.size1(), M.size2());
//	boost::numeric::ublas::vector<double> e(M.size1());
//	double duration = 0;
//	//NAME:  jacobiPseudoAsync
//	writeToAllStreams("\nJacobiAsyncTest", fp_outs);
//	// BEGIN TEST
//	begin = high_resolution_clock::now();
//	jacobiAsync(M, e, U, iter);
//	end = high_resolution_clock::now();
//	// END TEST
//
//	duration = duration_cast<nanoseconds>(end - begin).count() / 1000000.0;
//
//	// INFO
//	if (isWriteToConsole == "true") {
//		writeToAllStreams((boost::format("Eigenvalues: %1% \n U: %2% \nIter %3%\n Elapsed: %4%")
//			% e %U%iter%duration).str(), fp_outs);
//	}
//	else {
//		writeToAllStreams("Elapsed " + std::to_string(duration), fp_outs);
//	}
//
//	writeToAllStreams("============================", fp_outs);
//}

//SSTEQR
void LapackTest(boost::numeric::ublas::matrix<double> M, std::string isWriteToConsole, std::ofstream fp_outs[1], int i) {
	int iter;
	auto begin = std::chrono::high_resolution_clock::now();
	auto end = std::chrono::high_resolution_clock::now();
	if (isWriteToConsole == "true") {
		writeToAllStreams((boost::format("#%1%: \n") % i).str(), fp_outs);
		//writeToAllStreams((boost::format("A%1%: \n %2%") % i %M).str(), fp_outs);
	}
	else {
		writeToAllStreams(std::to_string(M.size1()) + "x" + std::to_string(M.size2()), fp_outs);
	}
	// PREP
	integer currMatrixSize = M.size1();
	real *diagonal = new real[currMatrixSize];
	real* offdiagonal = new real[currMatrixSize];
	
	//real diagonal[16], offdiagonal[15];

	for (int j = 0; j < currMatrixSize; j++) {
		diagonal[j] = M(j, j);
	}
	for (int j = 0; j < currMatrixSize - 1; j++) {
		offdiagonal[j] = M(j + 1, j);
	}
	real* dummy = new real[currMatrixSize*2];
	integer info = 0;
	// BEGIN TEST
	writeToAllStreams("\nLapackTest", fp_outs);
	begin = std::chrono::high_resolution_clock::now();
	char c4;
	c4 = 'N';
	//ssteqr_("N", &currMatrixSize, diagonal, offdiagonal, dummy, &currMatrixSize, dummy, info);
	ssteqr_(&c4, &currMatrixSize, diagonal, offdiagonal, dummy, &currMatrixSize, dummy, &info);
	end = std::chrono::high_resolution_clock::now();
	// END TEST

	double duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1000000.0;

	// INFO
	if (isWriteToConsole == "true") {
		std::string eig = "[";
		eig+=std::to_string(currMatrixSize);
		eig+="](";
		for (int i = 0; i < currMatrixSize -1; i++)
		{
			eig += std::to_string(diagonal[i]);
			eig += ", ";
		}
		eig += std::to_string(diagonal[currMatrixSize-1]);
		eig += ")";
		writeToAllStreams((boost::format("Eigenvalues: %1% \n Elapsed: %2%")
			%eig%duration).str(), fp_outs);
	}
	else {

		writeToAllStreams("Elapsed " + std::to_string(duration), fp_outs);
	}

	writeToAllStreams("============================", fp_outs);

}
int main(int argc, char **argv)
{

	//TBD: review code for memory overhead
	int numberOfMatrix = 3;
	std::string isWriteToConsole = "true";
	if (argc > 1 && argv) {

		numberOfMatrix = std::stoi(argv[1]);
		isWriteToConsole = argv[2];
	}

	std::ofstream fp_outs[1];
	fp_outs[0].open("output.txt", std::ios::out);
	
	boost::numeric::ublas::matrix<double>*MatrixArray = readFromSample(numberOfMatrix, "input.txt");
	std::cout << "read completed" << std::endl;
	int shift = 2;
	for (int i = shift; i < shift + 1; i++)
	{
		//doesn't work 
		JacobiGivensSyncTest(MatrixArray[i], isWriteToConsole, fp_outs, i);
		LapackTest(MatrixArray[i], isWriteToConsole, fp_outs, i);
		//JacobiSyncTest(MatrixArray[i], isWriteToConsole, fp_outs, i);
		//doesn't work 
		//JacobiAsyncTest(MatrixArray[i], isWriteToConsole, fp_outs, i);

		//writeToAllStreams("\n\n--------------------------------", fp_outs);
	}

	std::cin.get();	
	//return exit_success;
	return 0;
	
}
