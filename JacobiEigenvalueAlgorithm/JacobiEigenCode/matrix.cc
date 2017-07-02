#include "matrix.h"

std::ostream& operator<<(std::ostream& os, const matrix& mat)
{
	os << "\n[";
	for (int i = 0; i<mat.size(); ++i)
	{
		os << "	";
		for (int j = 0; j<mat.size(); ++j)
		{
			if (mat(i, j) >= 0.0f) os << " ";
			os << std::setw(10)<< std::fixed << std::setprecision(4) << mat(i, j);
			if (j == mat.size())
				os << " ";
			else
				os << ", ";
		}
		os << ";\n";
	}
	os << "]\n";
	return os;
}

