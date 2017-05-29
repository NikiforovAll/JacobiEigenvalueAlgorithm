#pragma once
#include <iostream>
#include <iomanip>
#include <cmath>

//typedef int ptr_size;
typedef long ptr_size;

template<typename T>
class matrix_ptrimpl_aligned
{
public:
	typedef T value_type;
	typedef T& reference;
	typedef const T& const_reference;

	explicit matrix_ptrimpl_aligned(int size)
		: n(size), datasize(size)
	{
		// If n is odd allocate an even sized matrix and pad with zeros.
		if (datasize % 2 == 1) datasize++;
		initialise();
		// Zero last row and column for the case when n is odd.
		for (int i = 0; i<datasize; ++i) d[datasize - 1][i] = 0.0f;
		for (int i = 0; i<datasize; ++i) d[i][datasize - 1] = 0.0f;
	}
	matrix_ptrimpl_aligned(const matrix_ptrimpl_aligned& m)
		: n(m.n), datasize(m.datasize)
	{
		initialise();
		// Copy data.
		for (int i = 0; i<datasize; ++i)
			for (int j = 0; j<datasize; ++j)
				d[i][j] = m.get(i, j);
	}
	~matrix_ptrimpl_aligned()
	{
		for (int i = 0; i<datasize; ++i)
			delete[] p[i];
		delete[] p;
		delete[] d;
	}

	void transpose()
	{
		char** oldp = p;	// Raw data array.
		T** oldd = d;	// Usable area (aligned).

		initialise();

		for (int i = 0; i<datasize; ++i)
			for (int j = 0; j<datasize; ++j)
				d[i][j] = oldd[j][i];

		for (int i = 0; i<datasize; ++i)
			delete[] oldp[i];
		delete[] oldp;
		delete[] oldd;
	}

	inline reference operator()(int i, int j)
	{
		return d[i][j];
	}

	inline reference get(int i, int j)
	{
		return d[i][j];
	}

	inline const_reference get(int i, int j) const
	{
		return d[i][j];
	}

	inline const_reference operator()(int i, int j) const
	{
		return d[i][j];
	}

	inline T* get_row(int r)
	{
		return d[r];
	}

	inline const T* get_row(int r) const
	{
		return d[r];
	}

	inline int size() const
	{
		return datasize;
	}

	inline int actual_size() const
	{
		return n;
	}

	void swap_rows(int i, int j)
	{
		T* t = d[i];
		d[i] = d[j];
		d[j] = t;
	}
private:
	void initialise()
	{
		// Allocate an array of Ts for each row.
		d = new T*[datasize];
		p = new char*[datasize];
		for (int i = 0; i<datasize; ++i)
		{
			// Allocate actual data on 16 byte boundaries and access through
			// d[i] array.
			p[i] = new char[datasize * sizeof(T) + 16];
			d[i] = (T*)(((ptr_size)p[i] + 16) & 0xfffffff0);
		}
	}

	const int n;
	int datasize;
	char** p;	// Raw data array.
	T** d;	// Usable area (aligned).
};
typedef matrix_ptrimpl_aligned<float> matrix;

inline double frobenius_norm(const matrix& mat)
{
	typedef matrix::value_type value_type;
	const int n = mat.size();
	double mag = 0.0;
	for (int i = 0; i<n; ++i)
	{
		const value_type *a = mat.get_row(i);
		for (int j = 0; j<n; ++j, ++a)
			mag += (*a)*(*a);
	}
	return std::sqrt(mag);
}

// Note: Using double precision arithmetic to avoid overflows
inline double off_diagonal_magnitude(const matrix& mat)
{
	typedef matrix::value_type value_type;
	const int n = mat.size();
	double mag = 0.0;
	for (int i = 0; i<n; ++i)
	{
		const value_type *a = mat.get_row(i);
		for (int j = 0; j<n; ++j, ++a)
		{
			if (i != j) mag += (*a)*(*a);
		}
	}
	return std::sqrt(mag);
}

std::ostream& operator<<(std::ostream& os, const matrix& M);

