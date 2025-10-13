#ifndef QUANTERNION_H_
#define QUANTERNION_H_

#include <math.h>

class Quanternion
{
	private:

	public:
		float s, x, y, z;
		Quanternion() :
			s(0), x(0), y(0), z(0)
		{

		}

		Quanternion(float scalar, float i, float j, float k)
		: s(scalar), x(i), y(j), z(k)
		{

		}

		float GetNorm();
		void Normalise();
		void Conjugate();
		void Inverse();
};

#endif /* QUANTERNION_H_ */
