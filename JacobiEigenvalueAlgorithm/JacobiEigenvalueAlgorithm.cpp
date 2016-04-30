#include <iostream>
#include <omp.h>

using namespace std;

int main(int argc, char **argv)
{
	int test(999);

	omp_set_num_threads(2);
#pragma omp parallel reduction(+:test)
	{
#pragma omp critical
		cout << "test = " << test << endl;
	}

	return EXIT_SUCCESS;
}