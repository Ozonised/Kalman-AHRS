#ifndef QUATERNION_H_
#define QUATERNION_H_

#include <math.h>

class Quaternion
{
	private:

	public:
		float s, x, y, z;
		Quaternion() :
			s(0), x(0), y(0), z(0)
		{

		}

		Quaternion(float scalar, float i, float j, float k)
		: s(scalar), x(i), y(j), z(k)
		{

		}

		float GetNorm();
		void Normalise();
		void Conjugate();
		void Inverse();
};

#endif /* QUATERNION_H_ */
