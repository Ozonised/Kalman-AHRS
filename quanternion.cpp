#include "quanternion.hpp"

float Quanternion::GetNorm()
{
	return sqrt(s * s + x * x + y * y + z * z);
}

void Quanternion::Normalise()
{
	double norm = GetNorm();
	s = s / norm;
	x = x / norm;
	y = y / norm;
	z = z / norm;
}

void Quanternion::Conjugate()
{
	x = -x;
	y = -y;
	z = -z;
}

void Quanternion::Inverse()
{
	double normSqr = GetNorm();
	normSqr *= normSqr;

	s /= normSqr;
	x /= normSqr;
	y /= normSqr;
	z /= normSqr;
}
