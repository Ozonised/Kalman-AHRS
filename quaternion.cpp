#include "quaternion.hpp"

float Quaternion::GetNorm()
{
	return sqrt(s * s + x * x + y * y + z * z);
}

void Quaternion::Normalise()
{
	float norm = GetNorm();
	if (norm)
	{
		s = s / norm;
		x = x / norm;
		y = y / norm;
		z = z / norm;
	}
}

void Quaternion::Conjugate()
{
	x = -x;
	y = -y;
	z = -z;
}

void Quaternion::Inverse()
{
	float normSqr = GetNorm();
	if (normSqr)
	{
		normSqr *= normSqr;

		s /= normSqr;
		x /= normSqr;
		y /= normSqr;
		z /= normSqr;
	}
}
