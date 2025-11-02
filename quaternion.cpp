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

/**
 * @brief Updates the quaternion object from Euler angles in Z-Y-X sequence
 *
 * @param[in] roll roll angle in radians
 * @param[in] pitch pitch angle in radians
 * @param[in] yaw yaw angle in radians
 *
 * @return None
 *
 * @note
 * 		- The Euler angles should be in 3-2-1 (Z-Y-X) sequence
 * 		- Refer: https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
 */
void Quaternion::EulerToQuat(float roll, float pitch, float yaw)
{
    float cr = cos(roll * 0.5f);
    float sr = sin(roll * 0.5f);
    float cp = cos(pitch * 0.5f);
    float sp = sin(pitch * 0.5f);
    float cy = cos(yaw * 0.5f);
    float sy = sin(yaw * 0.5f);

    s = cr * cp * cy + sr * sp * sy;
    x = sr * cp * cy - cr * sp * sy;
    y = cr * sp * cy + sr * cp * sy;
    z = cr * cp * sy - sr * sp * cy;
}

/**
 * @brief Retrieves the Euler angles in Z-Y-X sequence from quaternion object
 *
 * @param[out] roll roll angle in radians
 * @param[out] pitch pitch angle in radians
 * @param[out] yaw yaw angle in radians
 *
 * @return None
 *
 * @note
 * 		- Use this after calling Run() to get valid results
 * 		- The Euler angles are in 3-2-1 (Z-Y-X) sequence
 * 		- Can Gimbal lock
 * 		- Refer: https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
 */
void Quaternion::QuatToEuler(float& roll, float& pitch, float& yaw)
{
    // roll (x-axis rotation)
    float sinr_cosp = 2.0f * (s * x + y * z);
    float cosr_cosp = 1.0f - 2.0f * (x * x + y * y);
    roll = atan2(sinr_cosp, cosr_cosp);

    // pitch (y-axis rotation)
    float sinp =  2.0f * (s * y - x * z);
    pitch = asin(sinp);

    // yaw (z-axis rotation)
    float siny_cosp = 2.0f * (s * z + x * y);
    float cosy_cosp = 1.0f - 2.0f * (y * y + z * z);
    yaw = atan2(siny_cosp, cosy_cosp);
}
