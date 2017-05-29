#pragma once

#include <vector>
#include <limits>
#include <cmath>
//#ifdef enable_openmp
#include <omp.h>
//#endif

#include "matrix.h"
//#include "timer.h"

namespace parallel_jacobi
{
	//==========================================================================
	// Computes the minimum and maximum of two elements of arbitrary types
	//==========================================================================
	template<typename T>
	inline void minmax(T a, T b, T& min, T& max) {
		if (a > b) {
			max = a;
			min = b;
		}
		else {
			min = a;
			max = b;
		}
	}

	//==========================================================================
	// Functions to determine convergence
	//==========================================================================
	class converge_off_threshold {
		const double t;
		double tlast;
	public:
		// eps = tol*||A||_F
		converge_off_threshold(const double tolerance, const matrix& mat)
			: t(tolerance * frobenius_norm(mat)) { }

		bool not_converged(matrix& mat) {
			tlast = off_diagonal_magnitude(mat);
			return tlast > t;
		}

		void print() {
			std::cout << "threshold " << tlast << "/" << t << "\n";
		}
	};

	class converge_max_iterations {
		int i;
		const int imax;
	public:
		converge_max_iterations(const int iterations)
			: imax(iterations), i(0) { }

		bool not_converged(matrix& mat) {
			return ++i < imax;
		}

		void print() {
			std::cout << "iteration " << i << "/" << imax << ".\n";
		}
	};

	class converge_off_difference {
		double lastnorm;
		double lastnorm2;
		const double t;
	public:
		converge_off_difference(const double tolerance)
			: t(tolerance)
			, lastnorm(std::numeric_limits<double>::max()) { }

		bool not_converged(matrix& mat) {
			double o = off_diagonal_magnitude(mat);

			if ((lastnorm - o)<t)
				return false;

			lastnorm2 = lastnorm;
			lastnorm = o;
			return true;
		}
		void print() {
			std::cout << "off-diagonal difference " << (lastnorm2 - lastnorm) << "/" << t << ".\n";
		}
	};

	class music_permutation {
		std::vector<int> top, bot;
		const int n;
	public:
		music_permutation(int n) : n(n) {
			const int m = n / 2;
			top.resize(m);
			bot.resize(m);
			for (int i = 0; i<m; ++i) {
				top[i] = 2 * i;
				bot[i] = 2 * i + 1;
			}
		}

		inline void get(int k, int& p, int& q) {
			minmax(top[k], bot[k], p, q);
		}

		void permute() {
			// no permutation possible for 2x2 matrix
			if (n == 2) return;
			int m = n / 2;
			// Store end element to move down after bottom row is shifted.
			int e = top[m - 1];
			// Cycle top elements right and bottom elements left, being careful
			// not to overwrite anthing important.
			top[0] = bot[0];
			for (int k = 0, l = m - 1; k<m - 1; ++k, --l) {
				bot[k] = bot[k + 1];
				top[l] = top[l - 1];
			}
			// Move stored end element down to bot.
			bot[m - 1] = e;
			// Fix first top value.
			top[0] = 0;
		}
	};

	//==========================================================================
	// Symmetric Schur decomposition
	// Source: Algorithm 8.4.1 from Golub/Van Loan p.428
	//==========================================================================
	template<typename T>
	void symmetric_schur_new(const matrix& A, const unsigned int p,
		const unsigned int q, T& c, T& s)
	{
		const T epsilon = 1e-6;
		if (fabs(A(p, q)) > epsilon)
		{
			T tau = (A(q, q) - A(p, p)) / (2.0 * A(p, q));
			T t;
			if (tau >= 0) t = 1.0 / (tau + sqrt(1.0 + tau*tau));
			else t = -1.0 / (-tau + sqrt(1.0 + tau*tau));
			c = 1.0 / sqrt(1.0 + t*t);
			s = t*c;
		}
		else
		{
			c = 1.0;
			s = 0.0;
		}
	}

	//==========================================================================
	// Pre-multiply Jacobi rotation matrix
	//==========================================================================
	template<typename T>
	inline void premultiply(matrix& mat, const int p, const int q, const T c,
		const T s)
	{
		typedef matrix::value_type value_type;
		const int n = mat.size();
		value_type *rowp = mat.get_row(p),
			*rowq = mat.get_row(q);
		for (int i = 0; i<n; ++i, ++rowp, ++rowq)
		{
			const value_type mpi = *rowp,
				mqi = *rowq;
			*rowp = c*mpi + -s*mqi;
			*rowq = s*mpi + c*mqi;
		}
	}

	//==========================================================================
	// Post-multiply Jacobi rotation matrix
	//==========================================================================
	template<typename T>
	inline void postmultiply(matrix& mat, const int p, const int q, const T c,
		const T s)
	{
		typedef matrix::value_type value_type;
		for (int i = 0; i<mat.size(); ++i)
		{
			const value_type mip = mat(i, p),
				miq = mat(i, q);
			mat(i, p) = c*mip + -s*miq;
			mat(i, q) = s*mip + c*miq;
		}
	}

	//==========================================================================
	// Run the parallel jacobi algorithm
	// mat - the matrix to operate on
	// sc - a class with a function not_converged(mat) to determine when to stop
	//  the algorithm.
	//==========================================================================
	template<class StoppingCriterion, class Permutation>
	void run(matrix& mat, StoppingCriterion& sc, Permutation& pe)
		//,timer& root)
	{
		typedef matrix::value_type value_type;

		const int n = mat.size();
		const int m = n / 2;

		// Store sine-cosine pairs as recalculating them for post-multiplication
		// gives different values causing the algorithm to fail.
		value_type* si = new value_type[m];
		value_type* co = new value_type[m];

		// Set up timers.
		/*timer* convergence = new timer("convergence");
		timer* permute = new timer("permute");
		timer* premult = new timer("pre-multiplication");
		timer* postmult = new timer("post-multiplication");
		root.add(convergence);
		root.add(premult);
		root.add(postmult);
		root.add(permute);
*/
		const bool isodd = (n > mat.actual_size());
		bool not_converged = true;

		//root.start();
		while (not_converged) {
			// Do n sets of non-conflicting rotations.
			for (int set = 0; set<n; ++set) {
				//std::cout << "set: " << set << "\n";
				//premult->start();
				// Pre-multiply mat with each of the non-conflicting Jacobi
				// rotation matrices in the set concurrently.
				//#ifdef enable_openmp
				#pragma omp parallel for default(none) shared(mat, n, si, co, pe, isodd)
				//#endif
				for (int k = 0; k<m; ++k) {
					int p, q;
					pe.get(k, p, q);

					// Skip this transformation if it is on the padding row of
					// an odd sized matrix.
					if (isodd && (p == n || q == n)) continue;

					symmetric_schur_new(mat, p, q, co[k], si[k]);
					premultiply(mat, p, q, co[k], si[k]);
				}
				//premult->end();

				//postmult->start();
				// Similarly, post-multiply mat with each of the non-conflicting
				// Jacobi rotation matrices in the set concurrently.
				//#ifdef enable_openmp
				#pragma omp parallel for default(none) shared(mat, n, si, co, pe, isodd)
				//#endif
				for (int k = 0; k<m; ++k) {
					int p, q;
					pe.get(k, p, q);

					// Skip this transformation if it is on the padding row of
					// an odd sized matrix.
					if (isodd && (p == n || q == n)) continue;

					postmultiply(mat, p, q, co[k], si[k]);
				}
				//postmult->end();

				// Calculate the next set of non-conflicting rotations.
				//permute->start();
				pe.permute();
				//permute->end();
			}

			//convergence->start();
			not_converged = sc.not_converged(mat);
			//convergence->end();

			// Print how near stopping condition is to being satisfied.
			//sc.print();
		}
		//root.end();

		delete[] si;
		delete[] co;
	}
}

