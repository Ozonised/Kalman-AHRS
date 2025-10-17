#include "ExtendedKalmanFilter.hpp"

ExtendedKalmanFilter::ExtendedKalmanFilter()
{
	q.s = 1;
	qcap.s = 1;

	memset(P, 0, sizeof(P));
	P[0][0] = 0.1;
	P[1][1] = 0.1;
	P[2][2] = 0.1;
	P[3][3] = 0.6;

	// default gyroscoppe noise
	SigmaOmega = 0.0025;

	// default accelerometer noise
	R[0] =0.01;
	R[1] =0.01;
	R[2] =0.01;

	// default magnetometer noise
	R[3] = 0.1;
	R[4] = 0.1;
	R[5] = 0.1;

	magDeclination = 0;
	rx = 1;
	rz = 0;
}

/**
 * @brief Sets the sampling time for the Extended Kalman Filter.
 *
 * @param[in] freq  Sampling frequency in Hertz (Hz).
 *
 * @return None
 *
 * @note
 * Ensure that `freq` is non-zero to avoid division by zero.
 */
void ExtendedKalmanFilter::SetSampleTime(float freq)
{
	dt = 1.0f / freq;
	dt2 = dt / 2;
}

/**
 * @brief Sets the magnetic dip angle
 *
 * @param[in] degrees Dip angle in degrees
 *
 * @return None
 */
void ExtendedKalmanFilter::SetMagneticDip(float degrees)
{
	float rad = degrees *  0.01745f;
	rx = cos(rad);
	rz = sin(rad);
}

/**
 * @brief Sets the gyroscope process noise level for the Extended Kalman Filter.
 *
 * @param[in] Noise  Gyroscope process noise variance (typically in (rad/s)²).
 *
 * This approximation assumes identical and isotropic noise across all
 * rotational axes and neglects cross-axis correlations. It greatly simplifies
 * computation while providing adequate performance for most embedded AHRS
 * and orientation-tracking applications.
 *
 * @ return None
 */
void ExtendedKalmanFilter::SetGyroNoise(float Noise)
{
    SigmaOmega = Noise;
}

/**
 * @brief Sets the measurement noise covariance values for the Extended Kalman Filter.
 *
 * This function defines the expected noise levels for both the accelerometer
 * and magnetometer measurements.
 *
 * @param[in] NoiseAx  Accelerometer noise variance along the X-axis.
 * @param[in] NoiseAy  Accelerometer noise variance along the Y-axis.
 * @param[in] NoiseAz  Accelerometer noise variance along the Z-axis.
 * @param[in] NoiseMx  Magnetometer noise variance along the X-axis.
 * @param[in] NoiseMy  Magnetometer noise variance along the Y-axis.
 * @param[in] NoiseMz  Magnetometer noise variance along the Z-axis.
 *
 * @note
 * - All noise values should be positive.
 * - Units typically correspond to (sensor units)², e.g., (m/s²)² for accelerometer,
 *   and (µT)² for magnetometer.
 */
void ExtendedKalmanFilter::SetR(float NoiseAx, float NoiseAy, float NoiseAz,
		float NoiseMx, float NoiseMy, float NoiseMz)
{
	R[0] = NoiseAx;
	R[1] = NoiseAy;
	R[2] = NoiseAz;
	R[3] = NoiseMx;
	R[4] = NoiseMy;
	R[5] = NoiseMz;
}

/**
 * @brief Executes one iteration of the Extended Kalman Filter to estimate orientation.
 *
 * This function performs both the prediction and update steps of the
 * Extended Kalman Filter (EKF) using the latest sensor readings from the
 * accelerometer, gyroscope, and magnetometer. It computes the estimated
 * orientation as a unit quaternion `q`.
 *
 * @param[in] ax  Accelerometer X-axis measurement (m/s²)
 * @param[in] ay  Accelerometer Y-axis measurement (m/s²)
 * @param[in] az  Accelerometer Z-axis measurement (m/s²)
 * @param[in] gx  Gyroscope X-axis measurement (dps/s)
 * @param[in] gy  Gyroscope Y-axis measurement (dps/s)
 * @param[in] gz  Gyroscope Z-axis measurement (dps/s)
 * @param[in] mx  Magnetometer X-axis measurement (µT)
 * @param[in] my  Magnetometer Y-axis measurement (µT)
 * @param[in] mz  Magnetometer Z-axis measurement (µT)
 *
 * @return `true` if the computation completes successfully (no NaN or invalid values);
 *         `false` if any numerical instability or invalid computation occurs.
 * @note
 * - The function assumes sensor axes are aligned and calibrated.
 * - Ensure a valid sampling time is set via `SetSampleTime()` before calling `Run()`.
 * - If NaN or invalid results occur, the function returns `false`.
 */
bool ExtendedKalmanFilter::Run(float ax, float ay, float az, float gx, float gy, float gz, float mx, float my, float mz)
{
	// Normalize accelerometer and magnetometer readings
	float normA = sqrtf(ax * ax + ay * ay + az * az);
	float normM = sqrtf(mx * mx + my * my + mz * mz);

	if (normA)
	{
		ax = ax / normA;	ay = ay / normA;	az = az / normA;
	}

	if (normM)
	{
		mx = mx / normM;	my = my / normM; mz = mz / normM;
	}

	// convert from degrees/s to rads/s
	gx *= 0.01745f;
	gy *= 0.01745f;
	gz *= 0.01745f;

    // ============================================================
    // PREDICTION STEP - Using gyroscope
    // ============================================================

	qcap.s = q.s - dt2 * (gx*q.x + gy*q.y + gz*q.z);
	qcap.x = q.x + dt2 * (gx*q.s + gz*q.y - gy*q.z);
	qcap.y = q.y + dt2 * (gy*q.s + gx*q.z - gz*q.x);
	qcap.z = q.z + dt2 * (gz*q.s + gy*q.x - gx*q.y);

	qcap.Normalise();

	// Fill F
    F[0][0] = 1.0;         F[0][1] = -dt2 * gx;  	F[0][2] = -dt2 * gy;  	F[0][3] = -dt2 * gz;
    F[1][0] = dt2 * gx;  F[1][1] = 1.0;         	F[1][2] = dt2 * gz;   	F[1][3] = -dt2 * gy;
    F[2][0] = dt2 * gy;  F[2][1] = -dt2 * gz; 	F[2][2] = 1.0;          	F[2][3] = dt2 * gx;
    F[3][0] = dt2 * gz;  F[3][1] = dt2 * gy;  	F[3][2] = -dt2 * gx;  	F[3][3] = 1.0;

	// Process Noise Covariance: Pcap = FPFt + Q
	// FP = F * P
    float FP[4][4] = {0};
	// Pcap = F * P
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			for (int k = 0; k < 4; k++)
			{
				if (isnan(P[k][j]))
					return false;
				FP[i][j] += F[i][k] * P[k][j];
			}

	// Finally, Pcap = F*P*Ft + Q, here Q = SigmaOmega, i.e, gyro noise
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
		{
			Pcap[i][j] = 0;
			for (int k = 0; k < 4; k++)
				Pcap[i][j] += FP[i][k] * F[j][k];
			if (isnan(P[i][j]))
				return false;
		}

	Pcap[0][0] += SigmaOmega;
	Pcap[1][1] += SigmaOmega;
	Pcap[2][2] += SigmaOmega;
	Pcap[3][3] += SigmaOmega;

    // ============================================================
    // MEASUREMENT UPDATE
    // ============================================================
	float qcapS_Squared = qcap.s * qcap.s;
	float qcapX_Squared = qcap.x * qcap.x;
	float qcapY_Squared = qcap.y * qcap.y;
	float qcapZ_Squared = qcap.z * qcap.z;

	// Predicted accelerometer measurement (gravity in body frame)
	// acap = Cqt * g (NED frame)
	acap[0] = 2*(qcap.x*qcap.z+qcap.s*qcap.y);
	acap[1] = 2*(qcap.y*qcap.z-qcap.s*qcap.x);
	acap[2] = qcapS_Squared-qcapX_Squared-qcapY_Squared+qcapZ_Squared;


	if (normM)
	{
		// mcap = Cqt * r (NED frame)
		mcap[0] = (qcapS_Squared+qcapX_Squared-qcapY_Squared-qcapZ_Squared)*rx+2.0f*(qcap.x*qcap.z-qcap.s*qcap.y)*rz;
		mcap[1] = 2.0f*(qcap.x*qcap.y+qcap.s*qcap.z)*rx+2.0f*(qcap.y*qcap.z-qcap.s*qcap.x)*rz;
		mcap[2] = 2.0f*(qcap.x*qcap.z+qcap.s*qcap.y)*rx+(qcapS_Squared-qcapX_Squared-qcapY_Squared+qcapZ_Squared)*rz;
	}
	// v = z - h(qcap)
	// z = [a, m]^T and h(qcap) = [acap, mcap]^T
	v[0] = ax - acap[0];
	v[1] = ay - acap[1];
	v[2] = az - acap[2];

	if (normM)
	{
		v[3] = mx - mcap[0];
		v[4] = my - mcap[1];
		v[5] = mz - mcap[2];
	} else
	{
        // If no magnetometer, zero out these residuals to ignore magnetometer update
        v[3] = 0.0f;
        v[4] = 0.0f;
        v[5] = 0.0f;
	}

	// S = H(qcap)*Pcap*H(qcap)^T + R
	// H(qcap) = Jacobian of h(qcap)
	Hqcap[0][0] = 2.0f * qcap.y;    Hqcap[0][1] = 2.0f * qcap.z;    Hqcap[0][2] = 2.0f * qcap.s;    Hqcap[0][3] = 2.0f * qcap.x;
	Hqcap[1][0] = -2.0f * qcap.x;   Hqcap[1][1] = -2.0f * qcap.s;   Hqcap[1][2] = 2.0f * qcap.z;    Hqcap[1][3] = 2.0f * qcap.y;
	Hqcap[2][0] = 2.0f * qcap.s;    Hqcap[2][1] = -2.0f * qcap.x;   Hqcap[2][2] = -2.0f * qcap.y;   Hqcap[2][3] = 2.0f * qcap.z;

	Hqcap[3][0] = 2.0f * (rx*qcap.s - rz*qcap.y);  Hqcap[3][1] = 2.0f * (rx*qcap.x + rz*qcap.z);  Hqcap[3][2] = 2.0f * (-rx*qcap.y - rz*qcap.s); Hqcap[3][3] = 2.0f * (rx*qcap.z - rz*qcap.x);
	Hqcap[4][0] = 2.0f * (rx*qcap.z + rz*qcap.x);  Hqcap[4][1] = 2.0f * (rx*qcap.y + rz*qcap.s);  Hqcap[4][2] = 2.0f * (rx*qcap.x - rz*qcap.y);  Hqcap[4][3] = 2.0f * (rx*qcap.s + rz*qcap.z);
	Hqcap[5][0] = 2.0f * (-rx*qcap.y + rz*qcap.s); Hqcap[5][1] = 2.0f * (rx*qcap.z - rz*qcap.x);  Hqcap[5][2] = 2.0f * (rx*qcap.s + rz*qcap.y);  Hqcap[5][3] = 2.0f * (rx*qcap.y - rz*qcap.s);

	// PcapHT = Pcap*H(qcap)^T
	PcapHT[0][0] = Pcap[0][0]*Hqcap[0][0]+Pcap[0][1]*Hqcap[0][1]+Pcap[0][2]*Hqcap[0][2]+Pcap[0][3]*Hqcap[0][3];
	PcapHT[0][1] = Pcap[0][0]*Hqcap[1][0]+Pcap[0][1]*Hqcap[1][1]+Pcap[0][2]*Hqcap[1][2]+Pcap[0][3]*Hqcap[1][3];
	PcapHT[0][2] = Pcap[0][0]*Hqcap[2][0]+Pcap[0][1]*Hqcap[2][1]+Pcap[0][2]*Hqcap[2][2]+Pcap[0][3]*Hqcap[2][3];
	PcapHT[0][3] = Pcap[0][0]*Hqcap[3][0]+Pcap[0][1]*Hqcap[3][1]+Pcap[0][2]*Hqcap[3][2]+Pcap[0][3]*Hqcap[3][3];
	PcapHT[0][4] = Pcap[0][0]*Hqcap[4][0]+Pcap[0][1]*Hqcap[4][1]+Pcap[0][2]*Hqcap[4][2]+Pcap[0][3]*Hqcap[4][3];
	PcapHT[0][5] = Pcap[0][0]*Hqcap[5][0]+Pcap[0][1]*Hqcap[5][1]+Pcap[0][2]*Hqcap[5][2]+Pcap[0][3]*Hqcap[5][3];
	PcapHT[1][0] = Pcap[1][0]*Hqcap[0][0]+Pcap[1][1]*Hqcap[0][1]+Pcap[1][2]*Hqcap[0][2]+Pcap[1][3]*Hqcap[0][3];
	PcapHT[1][1] = Pcap[1][0]*Hqcap[1][0]+Pcap[1][1]*Hqcap[1][1]+Pcap[1][2]*Hqcap[1][2]+Pcap[1][3]*Hqcap[1][3];
	PcapHT[1][2] = Pcap[1][0]*Hqcap[2][0]+Pcap[1][1]*Hqcap[2][1]+Pcap[1][2]*Hqcap[2][2]+Pcap[1][3]*Hqcap[2][3];
	PcapHT[1][3] = Pcap[1][0]*Hqcap[3][0]+Pcap[1][1]*Hqcap[3][1]+Pcap[1][2]*Hqcap[3][2]+Pcap[1][3]*Hqcap[3][3];
	PcapHT[1][4] = Pcap[1][0]*Hqcap[4][0]+Pcap[1][1]*Hqcap[4][1]+Pcap[1][2]*Hqcap[4][2]+Pcap[1][3]*Hqcap[4][3];
	PcapHT[1][5] = Pcap[1][0]*Hqcap[5][0]+Pcap[1][1]*Hqcap[5][1]+Pcap[1][2]*Hqcap[5][2]+Pcap[1][3]*Hqcap[5][3];
	PcapHT[2][0] = Pcap[2][0]*Hqcap[0][0]+Pcap[2][1]*Hqcap[0][1]+Pcap[2][2]*Hqcap[0][2]+Pcap[2][3]*Hqcap[0][3];
	PcapHT[2][1] = Pcap[2][0]*Hqcap[1][0]+Pcap[2][1]*Hqcap[1][1]+Pcap[2][2]*Hqcap[1][2]+Pcap[2][3]*Hqcap[1][3];
	PcapHT[2][2] = Pcap[2][0]*Hqcap[2][0]+Pcap[2][1]*Hqcap[2][1]+Pcap[2][2]*Hqcap[2][2]+Pcap[2][3]*Hqcap[2][3];
	PcapHT[2][3] = Pcap[2][0]*Hqcap[3][0]+Pcap[2][1]*Hqcap[3][1]+Pcap[2][2]*Hqcap[3][2]+Pcap[2][3]*Hqcap[3][3];
	PcapHT[2][4] = Pcap[2][0]*Hqcap[4][0]+Pcap[2][1]*Hqcap[4][1]+Pcap[2][2]*Hqcap[4][2]+Pcap[2][3]*Hqcap[4][3];
	PcapHT[2][5] = Pcap[2][0]*Hqcap[5][0]+Pcap[2][1]*Hqcap[5][1]+Pcap[2][2]*Hqcap[5][2]+Pcap[2][3]*Hqcap[5][3];
	PcapHT[3][0] = Pcap[3][0]*Hqcap[0][0]+Pcap[3][1]*Hqcap[0][1]+Pcap[3][2]*Hqcap[0][2]+Pcap[3][3]*Hqcap[0][3];
	PcapHT[3][1] = Pcap[3][0]*Hqcap[1][0]+Pcap[3][1]*Hqcap[1][1]+Pcap[3][2]*Hqcap[1][2]+Pcap[3][3]*Hqcap[1][3];
	PcapHT[3][2] = Pcap[3][0]*Hqcap[2][0]+Pcap[3][1]*Hqcap[2][1]+Pcap[3][2]*Hqcap[2][2]+Pcap[3][3]*Hqcap[2][3];
	PcapHT[3][3] = Pcap[3][0]*Hqcap[3][0]+Pcap[3][1]*Hqcap[3][1]+Pcap[3][2]*Hqcap[3][2]+Pcap[3][3]*Hqcap[3][3];
	PcapHT[3][4] = Pcap[3][0]*Hqcap[4][0]+Pcap[3][1]*Hqcap[4][1]+Pcap[3][2]*Hqcap[4][2]+Pcap[3][3]*Hqcap[4][3];
	PcapHT[3][5] = Pcap[3][0]*Hqcap[5][0]+Pcap[3][1]*Hqcap[5][1]+Pcap[3][2]*Hqcap[5][2]+Pcap[3][3]*Hqcap[5][3];

	// S = H(qcap)*Pcap*H(qcap)^T + R
	// S = H(qcap)*PcapHT
	S[0][0] = Hqcap[0][0]*PcapHT[0][0]+Hqcap[0][1]*PcapHT[1][0]+Hqcap[0][2]*PcapHT[2][0]+Hqcap[0][3]*PcapHT[3][0];
	S[0][1] = Hqcap[0][0]*PcapHT[0][1]+Hqcap[0][1]*PcapHT[1][1]+Hqcap[0][2]*PcapHT[2][1]+Hqcap[0][3]*PcapHT[3][1];
	S[0][2] = Hqcap[0][0]*PcapHT[0][2]+Hqcap[0][1]*PcapHT[1][2]+Hqcap[0][2]*PcapHT[2][2]+Hqcap[0][3]*PcapHT[3][2];
	S[0][3] = Hqcap[0][0]*PcapHT[0][3]+Hqcap[0][1]*PcapHT[1][3]+Hqcap[0][2]*PcapHT[2][3]+Hqcap[0][3]*PcapHT[3][3];
	S[0][4] = Hqcap[0][0]*PcapHT[0][4]+Hqcap[0][1]*PcapHT[1][4]+Hqcap[0][2]*PcapHT[2][4]+Hqcap[0][3]*PcapHT[3][4];
	S[0][5] = Hqcap[0][0]*PcapHT[0][5]+Hqcap[0][1]*PcapHT[1][5]+Hqcap[0][2]*PcapHT[2][5]+Hqcap[0][3]*PcapHT[3][5];
	S[1][0] = Hqcap[1][0]*PcapHT[0][0]+Hqcap[1][1]*PcapHT[1][0]+Hqcap[1][2]*PcapHT[2][0]+Hqcap[1][3]*PcapHT[3][0];
	S[1][1] = Hqcap[1][0]*PcapHT[0][1]+Hqcap[1][1]*PcapHT[1][1]+Hqcap[1][2]*PcapHT[2][1]+Hqcap[1][3]*PcapHT[3][1];
	S[1][2] = Hqcap[1][0]*PcapHT[0][2]+Hqcap[1][1]*PcapHT[1][2]+Hqcap[1][2]*PcapHT[2][2]+Hqcap[1][3]*PcapHT[3][2];
	S[1][3] = Hqcap[1][0]*PcapHT[0][3]+Hqcap[1][1]*PcapHT[1][3]+Hqcap[1][2]*PcapHT[2][3]+Hqcap[1][3]*PcapHT[3][3];
	S[1][4] = Hqcap[1][0]*PcapHT[0][4]+Hqcap[1][1]*PcapHT[1][4]+Hqcap[1][2]*PcapHT[2][4]+Hqcap[1][3]*PcapHT[3][4];
	S[1][5] = Hqcap[1][0]*PcapHT[0][5]+Hqcap[1][1]*PcapHT[1][5]+Hqcap[1][2]*PcapHT[2][5]+Hqcap[1][3]*PcapHT[3][5];
	S[2][0] = Hqcap[2][0]*PcapHT[0][0]+Hqcap[2][1]*PcapHT[1][0]+Hqcap[2][2]*PcapHT[2][0]+Hqcap[2][3]*PcapHT[3][0];
	S[2][1] = Hqcap[2][0]*PcapHT[0][1]+Hqcap[2][1]*PcapHT[1][1]+Hqcap[2][2]*PcapHT[2][1]+Hqcap[2][3]*PcapHT[3][1];
	S[2][2] = Hqcap[2][0]*PcapHT[0][2]+Hqcap[2][1]*PcapHT[1][2]+Hqcap[2][2]*PcapHT[2][2]+Hqcap[2][3]*PcapHT[3][2];
	S[2][3] = Hqcap[2][0]*PcapHT[0][3]+Hqcap[2][1]*PcapHT[1][3]+Hqcap[2][2]*PcapHT[2][3]+Hqcap[2][3]*PcapHT[3][3];
	S[2][4] = Hqcap[2][0]*PcapHT[0][4]+Hqcap[2][1]*PcapHT[1][4]+Hqcap[2][2]*PcapHT[2][4]+Hqcap[2][3]*PcapHT[3][4];
	S[2][5] = Hqcap[2][0]*PcapHT[0][5]+Hqcap[2][1]*PcapHT[1][5]+Hqcap[2][2]*PcapHT[2][5]+Hqcap[2][3]*PcapHT[3][5];
	S[3][0] = Hqcap[3][0]*PcapHT[0][0]+Hqcap[3][1]*PcapHT[1][0]+Hqcap[3][2]*PcapHT[2][0]+Hqcap[3][3]*PcapHT[3][0];
	S[3][1] = Hqcap[3][0]*PcapHT[0][1]+Hqcap[3][1]*PcapHT[1][1]+Hqcap[3][2]*PcapHT[2][1]+Hqcap[3][3]*PcapHT[3][1];
	S[3][2] = Hqcap[3][0]*PcapHT[0][2]+Hqcap[3][1]*PcapHT[1][2]+Hqcap[3][2]*PcapHT[2][2]+Hqcap[3][3]*PcapHT[3][2];
	S[3][3] = Hqcap[3][0]*PcapHT[0][3]+Hqcap[3][1]*PcapHT[1][3]+Hqcap[3][2]*PcapHT[2][3]+Hqcap[3][3]*PcapHT[3][3];
	S[3][4] = Hqcap[3][0]*PcapHT[0][4]+Hqcap[3][1]*PcapHT[1][4]+Hqcap[3][2]*PcapHT[2][4]+Hqcap[3][3]*PcapHT[3][4];
	S[3][5] = Hqcap[3][0]*PcapHT[0][5]+Hqcap[3][1]*PcapHT[1][5]+Hqcap[3][2]*PcapHT[2][5]+Hqcap[3][3]*PcapHT[3][5];
	S[4][0] = Hqcap[4][0]*PcapHT[0][0]+Hqcap[4][1]*PcapHT[1][0]+Hqcap[4][2]*PcapHT[2][0]+Hqcap[4][3]*PcapHT[3][0];
	S[4][1] = Hqcap[4][0]*PcapHT[0][1]+Hqcap[4][1]*PcapHT[1][1]+Hqcap[4][2]*PcapHT[2][1]+Hqcap[4][3]*PcapHT[3][1];
	S[4][2] = Hqcap[4][0]*PcapHT[0][2]+Hqcap[4][1]*PcapHT[1][2]+Hqcap[4][2]*PcapHT[2][2]+Hqcap[4][3]*PcapHT[3][2];
	S[4][3] = Hqcap[4][0]*PcapHT[0][3]+Hqcap[4][1]*PcapHT[1][3]+Hqcap[4][2]*PcapHT[2][3]+Hqcap[4][3]*PcapHT[3][3];
	S[4][4] = Hqcap[4][0]*PcapHT[0][4]+Hqcap[4][1]*PcapHT[1][4]+Hqcap[4][2]*PcapHT[2][4]+Hqcap[4][3]*PcapHT[3][4];
	S[4][5] = Hqcap[4][0]*PcapHT[0][5]+Hqcap[4][1]*PcapHT[1][5]+Hqcap[4][2]*PcapHT[2][5]+Hqcap[4][3]*PcapHT[3][5];
	S[5][0] = Hqcap[5][0]*PcapHT[0][0]+Hqcap[5][1]*PcapHT[1][0]+Hqcap[5][2]*PcapHT[2][0]+Hqcap[5][3]*PcapHT[3][0];
	S[5][1] = Hqcap[5][0]*PcapHT[0][1]+Hqcap[5][1]*PcapHT[1][1]+Hqcap[5][2]*PcapHT[2][1]+Hqcap[5][3]*PcapHT[3][1];
	S[5][2] = Hqcap[5][0]*PcapHT[0][2]+Hqcap[5][1]*PcapHT[1][2]+Hqcap[5][2]*PcapHT[2][2]+Hqcap[5][3]*PcapHT[3][2];
	S[5][3] = Hqcap[5][0]*PcapHT[0][3]+Hqcap[5][1]*PcapHT[1][3]+Hqcap[5][2]*PcapHT[2][3]+Hqcap[5][3]*PcapHT[3][3];
	S[5][4] = Hqcap[5][0]*PcapHT[0][4]+Hqcap[5][1]*PcapHT[1][4]+Hqcap[5][2]*PcapHT[2][4]+Hqcap[5][3]*PcapHT[3][4];
	S[5][5] = Hqcap[5][0]*PcapHT[0][5]+Hqcap[5][1]*PcapHT[1][5]+Hqcap[5][2]*PcapHT[2][5]+Hqcap[5][3]*PcapHT[3][5];


	// H(qcap)*Pcap*H(qcap)^T + R = S + R
	// R is measurement noise of the accelerometer and magnetometer, it is a 6x6 diagonal matrix
	S[0][0] += R[0];
	S[1][1] += R[1];
	S[2][2] += R[2];
	S[3][3] += R[3];
	S[4][4] += R[4];
	S[5][5] += R[5];

    // ============================================================
    // MATRIX INVERSION - Compute S^-1
    // ============================================================

	// calculating inverse of S matrix, using LU decomposition
	// S * S^-1 = I, using LU decomposition S = LU
	memcpy(s, S, sizeof(S));
	for (uint8_t p = 0; p < 5; ++p)
	{
		// singular pivot
		if (fabs(s[p][p]) <= 0.00001f)
			return false;

		for (uint8_t ur = p + 1; ur < 6; ++ur)
		{
			//updating the l matrix :l[ur][p]
			s[ur][p] /= s[p][p];

			for (uint8_t uc = p + 1; uc < 6; ++uc)
			{
				s[ur][uc] -= s[ur][p] * s[p][uc];
			}
		}
	}

	// L * U * S^-1 = I, now let U * S^-1 = Y
	// therefore, L * Y = I
	float Icol[6], Y[6];
 	for (uint8_t col = 0; col < 6; col++)
	{
 		// now, L * Y = I, solve for Y, where L, Y and I are an nxn matrix (n = 6)
 		memset(Icol, 0, sizeof(Icol));
 		memset(Y, 0, sizeof(Y));
 		Icol[col] = 1;

 		float sum = 0;
 		for (uint8_t i = 0; i < 6; i++)
 		{
 			sum = 0;
 			for (uint8_t c = 0; c < i; c++)
 				sum += s[i][c] * Y[c];
 			Y[i] = Icol[i] - sum;
 		}

 		// Since, U * X = Y, here X = one column of S^-1
        // Now, solve for S^-1 (sinv)
 		for (int8_t i = 5; i > -1; i--)
 		{
 			S[i][col] = Y[i];
 			for (int8_t j = i + 1; j < 6; j++)
 				S[i][col] -= s[i][j] * S[j][col];

 			// singular pivot
 			if (fabs(s[i][i]) <= 0.00001f)
 				return false;

 			S[i][col] /= s[i][i];
 		}
	}

    // ============================================================
     // KALMAN GAIN: Kg = Pcap * H^T * S^-1 (4x3 matrix)
     // ============================================================
 	// Kg = Pcap * HT * S^-1(sinv)
 	// here, Kg = PcapHT * sinv
 	for (uint8_t m = 0; m < 4; m++)
 		for (uint8_t p = 0; p < 6; p++)
 		{
 			Kg[m][p] = 0;
 			for (uint8_t n = 0; n < 6; n++)
 				Kg[m][p] += PcapHT[m][n] * S[n][p];
			if (isnan(Kg[m][p]))
				return false;
 		}

    // ============================================================
    // STATE CORRECTION: q = qcap + Kg * v
    // ============================================================

	// Correction: q = qcap + Kg(z - h(q))
 	// q = qcap + Kg * V
	q.s = qcap.s + (Kg[0][0]*v[0]+Kg[0][1]*v[1]+Kg[0][2]*v[2]+Kg[0][3]*v[3]+Kg[0][4]*v[4]+Kg[0][5]*v[5]);
	q.x = qcap.x + (Kg[1][0]*v[0]+Kg[1][1]*v[1]+Kg[1][2]*v[2]+Kg[1][3]*v[3]+Kg[1][4]*v[4]+Kg[1][5]*v[5]);
	q.y = qcap.y + (Kg[2][0]*v[0]+Kg[2][1]*v[1]+Kg[2][2]*v[2]+Kg[2][3]*v[3]+Kg[2][4]*v[4]+Kg[2][5]*v[5]);

	if (normM)
		q.z = qcap.z + (Kg[3][0]*v[0]+Kg[3][1]*v[1]+Kg[3][2]*v[2]+Kg[3][3]*v[3]+Kg[3][4]*v[4]+Kg[3][5]*v[5]);

    // ============================================================
    // COVARIANCE UPDATE: P = (I - Kg * H) * Pcap
    // ============================================================

 	// P = (I4 - Kg * H(qcap)) * Pcap
 	// Solving I_KH = (I4 - Kg * Hqcap)
 	float I_KH[4][4];
 	I_KH[0][0] = 1.0 - (Kg[0][0]*Hqcap[0][0]+Kg[0][1]*Hqcap[1][0]+Kg[0][2]*Hqcap[2][0]+Kg[0][3]*Hqcap[3][0]+Kg[0][4]*Hqcap[4][0]+Kg[0][5]*Hqcap[5][0]);
 	I_KH[0][1] = -(Kg[0][0]*Hqcap[0][1]+Kg[0][1]*Hqcap[1][1]+Kg[0][2]*Hqcap[2][1]+Kg[0][3]*Hqcap[3][1]+Kg[0][4]*Hqcap[4][1]+Kg[0][5]*Hqcap[5][1]);
 	I_KH[0][2] = -(Kg[0][0]*Hqcap[0][2]+Kg[0][1]*Hqcap[1][2]+Kg[0][2]*Hqcap[2][2]+Kg[0][3]*Hqcap[3][2]+Kg[0][4]*Hqcap[4][2]+Kg[0][5]*Hqcap[5][2]);
 	I_KH[0][3] = -(Kg[0][0]*Hqcap[0][3]+Kg[0][1]*Hqcap[1][3]+Kg[0][2]*Hqcap[2][3]+Kg[0][3]*Hqcap[3][3]+Kg[0][4]*Hqcap[4][3]+Kg[0][5]*Hqcap[5][3]);
 	I_KH[1][0] = -(Kg[1][0]*Hqcap[0][0]+Kg[1][1]*Hqcap[1][0]+Kg[1][2]*Hqcap[2][0]+Kg[1][3]*Hqcap[3][0]+Kg[1][4]*Hqcap[4][0]+Kg[1][5]*Hqcap[5][0]);
 	I_KH[1][1] = 1.0 - (Kg[1][0]*Hqcap[0][1]+Kg[1][1]*Hqcap[1][1]+Kg[1][2]*Hqcap[2][1]+Kg[1][3]*Hqcap[3][1]+Kg[1][4]*Hqcap[4][1]+Kg[1][5]*Hqcap[5][1]);
 	I_KH[1][2] = -(Kg[1][0]*Hqcap[0][2]+Kg[1][1]*Hqcap[1][2]+Kg[1][2]*Hqcap[2][2]+Kg[1][3]*Hqcap[3][2]+Kg[1][4]*Hqcap[4][2]+Kg[1][5]*Hqcap[5][2]);
 	I_KH[1][3] = -(Kg[1][0]*Hqcap[0][3]+Kg[1][1]*Hqcap[1][3]+Kg[1][2]*Hqcap[2][3]+Kg[1][3]*Hqcap[3][3]+Kg[1][4]*Hqcap[4][3]+Kg[1][5]*Hqcap[5][3]);
 	I_KH[2][0] = -(Kg[2][0]*Hqcap[0][0]+Kg[2][1]*Hqcap[1][0]+Kg[2][2]*Hqcap[2][0]+Kg[2][3]*Hqcap[3][0]+Kg[2][4]*Hqcap[4][0]+Kg[2][5]*Hqcap[5][0]);
 	I_KH[2][1] = -(Kg[2][0]*Hqcap[0][1]+Kg[2][1]*Hqcap[1][1]+Kg[2][2]*Hqcap[2][1]+Kg[2][3]*Hqcap[3][1]+Kg[2][4]*Hqcap[4][1]+Kg[2][5]*Hqcap[5][1]);
 	I_KH[2][2] = 1.0 - (Kg[2][0]*Hqcap[0][2]+Kg[2][1]*Hqcap[1][2]+Kg[2][2]*Hqcap[2][2]+Kg[2][3]*Hqcap[3][2]+Kg[2][4]*Hqcap[4][2]+Kg[2][5]*Hqcap[5][2]);
 	I_KH[2][3] = -(Kg[2][0]*Hqcap[0][3]+Kg[2][1]*Hqcap[1][3]+Kg[2][2]*Hqcap[2][3]+Kg[2][3]*Hqcap[3][3]+Kg[2][4]*Hqcap[4][3]+Kg[2][5]*Hqcap[5][3]);
 	I_KH[3][0] = -(Kg[3][0]*Hqcap[0][0]+Kg[3][1]*Hqcap[1][0]+Kg[3][2]*Hqcap[2][0]+Kg[3][3]*Hqcap[3][0]+Kg[3][4]*Hqcap[4][0]+Kg[3][5]*Hqcap[5][0]);
 	I_KH[3][1] = -(Kg[3][0]*Hqcap[0][1]+Kg[3][1]*Hqcap[1][1]+Kg[3][2]*Hqcap[2][1]+Kg[3][3]*Hqcap[3][1]+Kg[3][4]*Hqcap[4][1]+Kg[3][5]*Hqcap[5][1]);
 	I_KH[3][2] = -(Kg[3][0]*Hqcap[0][2]+Kg[3][1]*Hqcap[1][2]+Kg[3][2]*Hqcap[2][2]+Kg[3][3]*Hqcap[3][2]+Kg[3][4]*Hqcap[4][2]+Kg[3][5]*Hqcap[5][2]);
 	I_KH[3][3] = 1.0 - (Kg[3][0]*Hqcap[0][3]+Kg[3][1]*Hqcap[1][3]+Kg[3][2]*Hqcap[2][3]+Kg[3][3]*Hqcap[3][3]+Kg[3][4]*Hqcap[4][3]+Kg[3][5]*Hqcap[5][3]);

 	// Now, P = I_KH * Pcap
 	P[0][0] = I_KH[0][0]*Pcap[0][0]+I_KH[0][1]*Pcap[1][0]+I_KH[0][2]*Pcap[2][0]+I_KH[0][3]*Pcap[3][0];
 	P[0][1] = I_KH[0][0]*Pcap[0][1]+I_KH[0][1]*Pcap[1][1]+I_KH[0][2]*Pcap[2][1]+I_KH[0][3]*Pcap[3][1];
 	P[0][2] = I_KH[0][0]*Pcap[0][2]+I_KH[0][1]*Pcap[1][2]+I_KH[0][2]*Pcap[2][2]+I_KH[0][3]*Pcap[3][2];
 	P[0][3] = I_KH[0][0]*Pcap[0][3]+I_KH[0][1]*Pcap[1][3]+I_KH[0][2]*Pcap[2][3]+I_KH[0][3]*Pcap[3][3];
 	P[1][0] = I_KH[1][0]*Pcap[0][0]+I_KH[1][1]*Pcap[1][0]+I_KH[1][2]*Pcap[2][0]+I_KH[1][3]*Pcap[3][0];
 	P[1][1] = I_KH[1][0]*Pcap[0][1]+I_KH[1][1]*Pcap[1][1]+I_KH[1][2]*Pcap[2][1]+I_KH[1][3]*Pcap[3][1];
 	P[1][2] = I_KH[1][0]*Pcap[0][2]+I_KH[1][1]*Pcap[1][2]+I_KH[1][2]*Pcap[2][2]+I_KH[1][3]*Pcap[3][2];
 	P[1][3] = I_KH[1][0]*Pcap[0][3]+I_KH[1][1]*Pcap[1][3]+I_KH[1][2]*Pcap[2][3]+I_KH[1][3]*Pcap[3][3];
 	P[2][0] = I_KH[2][0]*Pcap[0][0]+I_KH[2][1]*Pcap[1][0]+I_KH[2][2]*Pcap[2][0]+I_KH[2][3]*Pcap[3][0];
 	P[2][1] = I_KH[2][0]*Pcap[0][1]+I_KH[2][1]*Pcap[1][1]+I_KH[2][2]*Pcap[2][1]+I_KH[2][3]*Pcap[3][1];
 	P[2][2] = I_KH[2][0]*Pcap[0][2]+I_KH[2][1]*Pcap[1][2]+I_KH[2][2]*Pcap[2][2]+I_KH[2][3]*Pcap[3][2];
 	P[2][3] = I_KH[2][0]*Pcap[0][3]+I_KH[2][1]*Pcap[1][3]+I_KH[2][2]*Pcap[2][3]+I_KH[2][3]*Pcap[3][3];
 	P[3][0] = I_KH[3][0]*Pcap[0][0]+I_KH[3][1]*Pcap[1][0]+I_KH[3][2]*Pcap[2][0]+I_KH[3][3]*Pcap[3][0];
 	P[3][1] = I_KH[3][0]*Pcap[0][1]+I_KH[3][1]*Pcap[1][1]+I_KH[3][2]*Pcap[2][1]+I_KH[3][3]*Pcap[3][1];
 	P[3][2] = I_KH[3][0]*Pcap[0][2]+I_KH[3][1]*Pcap[1][2]+I_KH[3][2]*Pcap[2][2]+I_KH[3][3]*Pcap[3][2];
 	P[3][3] = I_KH[3][0]*Pcap[0][3]+I_KH[3][1]*Pcap[1][3]+I_KH[3][2]*Pcap[2][3]+I_KH[3][3]*Pcap[3][3];

 	q.Normalise();
 	return true;
}

/**
 * @brief Retrieves the current estimated orientation quaternion.
 *
 * @param[out] qState Reference to `Quanternion` object that will receive
 *        	  the current orientation (s, x, y, z components).
 *
 * @details
 * The quaternion `q` represents the current estimated orientation of
 * the system with respect to the global (Earth) reference frame.
 *
 * @note
 * - Ensure that the filter has been updated (via `Run()`) before calling
 *   this function to obtain valid orientation data.
 */
void ExtendedKalmanFilter::GetOrientation(Quanternion& qState)
{
	qState.s = q.s;
	qState.x = q.x;
	qState.y = q.y;
	qState.z = q.z;
}
