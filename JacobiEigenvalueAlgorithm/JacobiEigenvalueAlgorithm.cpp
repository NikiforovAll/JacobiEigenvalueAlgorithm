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

using namespace std;
using namespace boost::numeric::ublas;
using namespace std::chrono;

void JacobiGivensSyncTest(matrix<double> M, string isWriteToConsole, ofstream fp_outs[1], int i) {
	int iter;
	auto begin = std::chrono::high_resolution_clock::now();
	auto end = std::chrono::high_resolution_clock::now();
	if (isWriteToConsole == "true") {
		writeToAllStreams((boost::format("A%1%: \n %2%") % i %M).str(), fp_outs);
	}
	else {
		writeToAllStreams(std::to_string(M.size1()) + "x" + std::to_string(M.size2()), fp_outs);
	}
	matrix<double> U(M.size1(), M.size2());
	boost::numeric::ublas::vector<double> e(M.size1());
	double duration = 0;
	//NAME:  jacobiPseudoAsync
	writeToAllStreams("\nJacobiGivensSyncTest", fp_outs);
	// BEGIN TEST
	begin = high_resolution_clock::now();
	jacobiPseudoAsync(M, e, U, iter);
	end = high_resolution_clock::now();
	// END TEST

	duration = duration_cast<nanoseconds>(end - begin).count() / 1000000.0;

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

void JacobiSyncTest(matrix<double> M, string isWriteToConsole, ofstream fp_outs[1], int i) {
	int iter;
	auto begin = std::chrono::high_resolution_clock::now();
	auto end = std::chrono::high_resolution_clock::now();
	if (isWriteToConsole == "true") {
		writeToAllStreams((boost::format("A%1%: \n %2%") % i %M).str(), fp_outs);
	}
	else {
		writeToAllStreams(std::to_string(M.size1()) + "x" + std::to_string(M.size2()), fp_outs);
	}
	matrix<double> U(M.size1(), M.size2());
	boost::numeric::ublas::vector<double> e(M.size1());
	double duration = 0;
	//NAME:  jacobiPseudoAsync
	writeToAllStreams("\nJacobiSyncTest", fp_outs);
	// BEGIN TEST
	begin = high_resolution_clock::now();
	jacobiSync(M, e, U, iter, false);
	end = high_resolution_clock::now();
	// END TEST

	duration = duration_cast<nanoseconds>(end - begin).count() / 1000000.0;

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

void JacobiAsyncTest(matrix<double> M, string isWriteToConsole, ofstream fp_outs[1], int i) {
	int iter;
	auto begin = std::chrono::high_resolution_clock::now();
	auto end = std::chrono::high_resolution_clock::now();
	if (isWriteToConsole == "true") {
		writeToAllStreams((boost::format("A%1%: \n %2%") % i %M).str(), fp_outs);
	}
	else {
		writeToAllStreams(std::to_string(M.size1()) + "x" + std::to_string(M.size2()), fp_outs);
	}
	matrix<double> U(M.size1(), M.size2());
	boost::numeric::ublas::vector<double> e(M.size1());
	double duration = 0;
	//NAME:  jacobiPseudoAsync
	writeToAllStreams("\nJacobiAsyncTest", fp_outs);
	// BEGIN TEST
	begin = high_resolution_clock::now();
	jacobiAsync(M, e, U, iter);
	end = high_resolution_clock::now();
	// END TEST

	duration = duration_cast<nanoseconds>(end - begin).count() / 1000000.0;

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


int main(int argc, char **argv)
{

	int numberOfMatrix = 3;
	string isWriteToConsole = "false";
	if (argc > 1 && argv) {

		numberOfMatrix = stoi(argv[1]);
		isWriteToConsole = argv[2];
	}


	ofstream fp_outs[1];
	fp_outs[0].open("output.txt", ios::out);

	matrix<double>*MatrixArray = readFromSample(numberOfMatrix, "input.txt");
	for (int i = 1; i < 2; i++)
	{
		JacobiGivensSyncTest(MatrixArray[i], isWriteToConsole, fp_outs, i);
		//JacobiSyncTest(MatrixArray[i], isWriteToConsole, fp_outs, i);
		//JacobiAsyncTest(MatrixArray[i], isWriteToConsole, fp_outs, i);

		writeToAllStreams("\n\n--------------------------------", fp_outs);
	}
	cin.get();	
	return EXIT_SUCCESS;
}