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
using namespace std;
using namespace boost::numeric::ublas;
using namespace boost;


#define _EPS 10e-03

matrix<double>* readFromSample(int num, string filename) {
	string line;
	ifstream myfile(filename);
	int count = 0;
	string matrixStr;
	matrix<double>* result;
	result = new (nothrow) matrix<double>[num];
	bool isMatrix = false;
	if (myfile.is_open())
	{
		while (getline(myfile, line))
		{
			if (count < num)
				if (line == "#") {
					bool sizeState = true;
					bool isReadState = false;
					string sizeStr;
					string rowStr = "";
					int size[2];
					int currIndex[2] = { 0, 0 };
					for (int i = 1; i < matrixStr.size() - 1; i++) {
						if (sizeState) {
							if (matrixStr[i] == ']') {
								std::istringstream ss(sizeStr);
								std::string token;
								int currSize = 0;
								while (std::getline(ss, token, ',')) {
									size[currSize] = stoi(token);
									currSize++;
								}
								sizeState = false;
								result[count] = matrix<double>(size[0], size[1]);
								i++;
							}
							else {
								sizeStr += matrixStr[i];
							}
						}
						else {
							if (matrixStr[i] == '(') {
								isReadState = true;
								currIndex[1] = 0;
								continue;
							}
							if (isReadState && matrixStr[i] != ')') {
								rowStr += matrixStr[i];
							}
							else if (matrixStr[i] == ')') {
								std::istringstream ss(rowStr);
								std::string token;
								int currSize = 0;
								while (std::getline(ss, token, ',')) {
									result[count](currIndex[0], currSize) = stod(token);
									currSize++;
								}
								currIndex[0]++;
								isReadState = false;
								rowStr = "";
								i++;
							}
						}
					}
					matrixStr = "";
					count++;
					if (count > num) {
						break;
					}
				}
				else {
					matrixStr += line;
				}

		}
		myfile.close();

	}
	else cout << "Unable to open file";
	return result;
}
// abso

void abs(matrix<double> &M)
{
	int n = M.size1();

	for (int k = 0; k < n; k++)
	{
		for (int i = 0; i < n; i++)
		{
			M(k, i) = abs(M(k, i));
		}
	}
}
double sumOffDiagonal(matrix<double> &S) {
	double sum = 0;
	for (int i = 0; i < S.size1(); i++) {
		for (int j = i + 1; j < S.size2() - 1; j++)
		{
			sum += abs(S(i, j));
		}
	}
	return sum;
}
void writeToAllStreams(string str, ofstream stream[1]) {
	cout << str << endl;
	stream[0] << str << endl;
}

// finding pivot (off diagonal element)
void findMax(matrix<double> &M, int &row, int &col)
{
	double m = M(0, 0);
	row = col = 0;
	int n = M.size1();

	for (int k = 0; k < n; k++)
	{
		for (int l = k; l < n; l++)
		{
			if (M(k, l) > m)
			{
				row = k;
				col = l;
				m = M(k, l);
			}
		}
	}
}

// rotation 
void rotateRowCol(matrix<double> &S, matrix<double> &U, int row, int col)
{
	double t, c, s, theta;
	int n = S.size1();
	if (row == col) return;
	theta = (S(row, row) - S(col, col)) / (2 * S(row, col));
	if (theta < _EPS)
	{
		t = 1 / (abs(theta) + sqrt(theta*theta + 1));
	}
	else {
		t = 1 / abs(2 * theta);
	}
	if (theta < 0) t = -t;
	c = 1 / sqrt(t*t + 1);
	s = c*t;
	matrix<double> R(n, n);
	R = identity_matrix<double>(n, n);
	R(row, row) = R(col, col) = c; R(row, col) = s; R(col, row) = -s;
	S = prod(S, trans(R));
	S = prod(R, S);
	U = prod(U, R);
}

//also possible to modify algorithm by changing multiplication to tridiagonal matrix form
int jacobiSync(matrix<double> &S, boost::numeric::ublas::vector<double> &e, matrix<double>  &U, int &iter)
{
	iter = 0;
	int col, row;
	bool iterating = true;
	int n = S.size1();
	if (S.size2() != n)
	{
		return -1;
	}
	matrix<double> M(n, n);
	e = zero_vector<double>(n);
	U = identity_matrix<double>(n, n);
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
		rotateRowCol(S, U, row, col);

		//cout<<row<<" row|col "<<col <<" sum: "<< sumOffDiagonal(S) << endl;
		//if (Smax < _EPS * norm_frobenius(S)) iterating = false;
		if (sumOffDiagonal(S) < _EPS) iterating = false;
	}

	for (int i = 0; i < n; i++) e(i) = S(i, i);

	return 0;
}


void generateDisJointPairs(boost::numeric::ublas::vector<double> &top, boost::numeric::ublas::vector<double> &bot) {
	int sizeTop = top.size();
	int sizeBot = bot.size();
	boost::numeric::ublas::vector<double> newTop(sizeTop);
	boost::numeric::ublas::vector<double> newBot(sizeBot);
	for (int i = 0; i < sizeTop; i++)
	{
		if (i == 0) {
			newTop(i) = bot(0);
		}
		else if (i > 0) {
			newTop(i) = top(i - 1);
		}
		if (i == (sizeTop - 1)) {
			newBot(i) = top(i);
		}
		else {
			newBot(i) = bot(i + 1);
		}
	}
	top = newTop;
	bot = newBot;
}

void generateStartDisJointPair(boost::numeric::ublas::vector<double> &top, boost::numeric::ublas::vector<double> &bot) {
	for (int i = 0, j = 0; (i + j) < top.size() * 2;) {
		if ((i + j) % 2 == 0) {
			bot(i) = (i + 1) * 2;
			i++;
		}
		else {
			top(j) = 1 + j * 2;
			j++;
		}
	}
}
void cycleVector(boost::numeric::ublas::vector<double> &v) {
	boost::numeric::ublas::vector<double> newV(v.size());
	//newV(0) = v(v.size - 1);
	for (int i = 1; i < v.size(); i++) {
		newV(i) = v(i - 1);
	}
	newV(0) = v(v.size()-1);
	v = newV;
}
bool checkConjunction(boost::numeric::ublas::vector<double> &v, int i) {
	int j = i - 1;
	bool result = (v(i)==0);
	if(i-1>=0){
		result = result && (v(i - 1) == 0);
	}
	if (i + 1 < v.size() - 1) {
		result = result && (v(i + 1) == 0);
	}
	return result;
}
void setConjunction(boost::numeric::ublas::vector<double> &v, int i) {
	v(i) = 1;
	if (i - 1 >= 0) {
		v(i - 1) = 1;
	}
	if (i + 1 < v.size() - 1) {
		v(i + 1) = 1;
	}
}
void RotateColRowGivens(matrix<double> &S, int row, int col, int len) {
	double x, y, t, s, c;
	int p, q;
	p = row;
	q = col;
	x = S(p, p - 1);
	y = S(q, p - 1);
	t = sqrt(x*x + y*y);
	s = -y / t;
	c = x / t;
	S(p, p - 1) = t;
	S(q, p - 1) = 0;
	// left multiplication 
	for (int i = 1; i <= len; i++) {
		x = S(p, p + i - 1);
		y = S(q, p + i - 1);
		S(p, p + i - 1) = c*x - s*y;
		S(q, p + i - 1) = s*x + c*y;
	}
	S(p - 1, p) = t;
	S(p - 1, q) = 0;
	// right multiplication 
	for (int i = 1; i <= len; i++) {
		x = S(p + i - 1, p);
		y = S(p + i - 1, q);
		S(p + i - 1, p) = c*x - s*y;
		S(p + i - 1, q) = s*x + c*y;
	}
	
}
int jacobiPseudoAsync(matrix<double> &S, boost::numeric::ublas::vector<double> &e, matrix<double>  &U, int &iter) {
	int n = S.size1();
	iter = 0;
	//Jacobi rotation
	for (int pp = 0; pp < n - 2; pp++) {
		int p = pp;
		int q = p + 1;
		while (abs(S(p, q)) > _EPS*sqrt(S(q, q)*S(q, q) + S(p, p)*(p, p))) {
			iter++;
			double tau = (S(q, q) - S(p, p)) / (2 * S(p, q));
			double c, s, t;
			if (tau == 0) {
				c = 1 / sqrt(2);
				s = c;
			}else {
				t = (tau >= 0 ? 1.0 : -1.0)/(abs(tau)+sqrt(1+tau*tau));
				c = 1.0 / sqrt(1 + t*t);
				s = t*c;
			}
			double x, y;
			// left multiplication 
			for (int i = 1; i <=3; i++) {
				x = S(p, p + i - 1);
				y = S(q, p + i - 1);
				S(p, p + i - 1) = c*x - s*y;
				S(q, p + i - 1) = s*x + c*y;
			}
			// right multiplication 
			for (int i = 1; i <= 3; i++) {
				x = S(p + i - 1, p);
				y = S(p + i - 1, q);
				S(p + i - 1,p) = c*x - s*y;
				S(p + i - 1,q) = s*x + c*y;
			}
			//matrix<double> U(S.size1(), S.size2());
			//rotateRowCol(S, U, q, p);
			/*cout << S(p, q)<<endl<<S(q,p)<<endl;
			cout << S << endl;
			cout << p <<" p|q "<< q << endl;*/
			// removing redundant accumulation via Givens rotation
			for (int p = 1; p < n - 2; p++) {
				q = p + 1;
				if (abs(S(q, p - 1)) > _EPS*sqrt(S(q, q)*S(q, q) + S(p, p)*(p, p))) {
					RotateColRowGivens(S, p, q, 3);
				}
			}
			p = n - 2;
			q = p+1;
			//last step
			if (abs(S(p, q)) > _EPS*sqrt(S(q, q)*S(q, q) + S(p, p)*(p, p))) {
				RotateColRowGivens(S, p, q, 2);
			}
			p = pp;
			q = p + 1;
			//cout <<"iter: "<<endl<<S;
			//cin.get();
		}
	}
	return 0;
}



int jacobiAsync(matrix<double> &S, boost::numeric::ublas::vector<double> &e, matrix<double>  &U, int &iter) {
	iter = 0;
	int col, row;
	bool iterating = true;
	const int n = S.size1();
	if (S.size2() != n)
	{
		return -1;
	}
	matrix<double> M(n, n);
	U = identity_matrix<double>(n, n);
	boost::numeric::ublas::vector<double> top(S.size1() / 2);
	boost::numeric::ublas::vector<double> bot(S.size2() / 2);
	generateStartDisJointPair(top, bot);
	while (iterating)
	{
		iter++;
		//TBD

		boost::numeric::ublas::vector<double> distr_status(n);
		matrix<double> toProcess(n, 2);
		int processPointer = 0;
		
		for (int j = 0; j < n - 1; j++)
		{
			distr_status.clear();
			toProcess.clear();
			processPointer = 0;
			#pragma omp parallel for shared(distr_status, top, bot, S, U, processPointer) private(row,col) //num_threads(2)
				//#pragma omp parallel for shared(top,bot) private(row,col)
				for (int i = 0; i < S.size1() / 2; i++)
				{
					//cout << "pair" << endl << top << endl << bot << endl;
					string num_thread_str = "[" +std::to_string(omp_get_thread_num())+ "]";
					row = std::max(top(i), bot(i)) - 1;
					col = std::min(top(i), bot(i)) - 1;
					if (checkConjunction(distr_status, row) && checkConjunction(distr_status, col)) {
						cout << "d_status" << num_thread_str << distr_status << endl;
						setConjunction(distr_status, row);
						setConjunction(distr_status, col);
						rotateRowCol(S, U, col, row);
						cout <<num_thread_str<<"operated:" << col << " -col-row- " << row << endl;
					}
					else {
						toProcess(processPointer, 1) = row;
						toProcess(processPointer, 0) = col;
						#pragma omp atomic
						processPointer++;
						cout << num_thread_str << "not operated:" << col << " -col-row- " << row << endl;
					}				

				}				
			
			//cout << distr_status << endl;
			//cout << "pair" << endl << top << endl << bot << endl;
			//cout << "toProcess" << toProcess <<endl;
			//cout << "barrier point -----------" << endl;
			for (int i = 0; (toProcess(i, 0) != 0 && (toProcess(i, 1) != 0)); i++) {
				rotateRowCol(S, U, toProcess(i,0), toProcess(i,1));
			}
			generateDisJointPairs(top, bot);
			//cin.get();
		}

		//cout << "pair" << endl << top << endl << bot << endl;
		//cin.get();
		cycleVector(bot);
		cout << S << endl;
		cout <<"sumOffDiagonal"<< sumOffDiagonal(S)<<endl;
		if (sumOffDiagonal(S) < _EPS) iterating = false;
	}

	for (int i = 0; i < n; i++) e(i) = S(i, i);

	return 0;
}

int main(int argc, char **argv)
{
	
	int numberOfMatrix = 2;
	if (argc > 1 && argv) {

		numberOfMatrix = stoi(argv[1]);
	}
	/*int test(999);

	omp_set_num_threads(2);
#pragma omp parallel reduction(+:test)
	{
#pragma omp critical
		std::cout << "test = " << test << std::endl;
	}
	getchar();*/
	auto begin = std::chrono::high_resolution_clock::now();
	auto end = std::chrono::high_resolution_clock::now();
	ofstream fp_outs[1];
	fp_outs[0].open("output.txt", ios::out);
	matrix<double>*MatrixArray = readFromSample(numberOfMatrix, "input.txt");
	for (int i = 0; i < numberOfMatrix; i++)
	{
		//boost::timer t; t.elapsed()

		writeToAllStreams("============================", fp_outs);
		begin = std::chrono::high_resolution_clock::now();
		int iter;
		matrix<double> M = MatrixArray[i];
		writeToAllStreams((boost::format("A%1%: \n %2%") % i %M).str(), fp_outs);
		writeToAllStreams("Sync version", fp_outs);
		matrix<double> U(M.size1(), M.size2());
		boost::numeric::ublas::vector<double> e(M.size1());

		jacobiSync(M, e, U, iter);
		end = std::chrono::high_resolution_clock::now();
		double duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1000000.0;
		writeToAllStreams((boost::format("Eigenvalues: %1% \n U: %2% \nIter %3%\n Elapsed: %4%")
			% e %U%iter%duration).str(), fp_outs);
		writeToAllStreams("============================", fp_outs);
		writeToAllStreams("Pseudo Sync version", fp_outs);

		M = MatrixArray[i];
		begin = std::chrono::high_resolution_clock::now();
		jacobiPseudoAsync(M, e, U, iter);
		end = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1000000.0;
		writeToAllStreams((boost::format("Eigenvalues: %1% \n U: %2% \nIter %3%\n Elapsed: %4%")
			% e %U%iter%duration).str(), fp_outs);
	/*	writeToAllStreams("Async version", fp_outs);
		M = MatrixArray[i];
		begin = std::chrono::high_resolution_clock::now();
		jacobiAsync(M, e, U, iter);
		end = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1000000.0;
		writeToAllStreams((boost::format("Eigenvalues: %1% \n U: %2% \nIter %3%\n Elapsed: %4%")
			% e %U%iter%duration).str(), fp_outs);*/
	}
	cin.get();
	//TBD: async solution
	return EXIT_SUCCESS;
}