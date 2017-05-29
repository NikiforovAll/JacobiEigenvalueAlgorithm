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
			os << std::fixed << std::setprecision(32) << mat(i, j);
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

