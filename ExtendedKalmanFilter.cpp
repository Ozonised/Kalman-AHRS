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
	SigmaOmega[0] = 0.0025f;
	SigmaOmega[1] = 0.0025f;
	SigmaOmega[2] = 0.0025f;

	// default accelerometer noise
	R[0] =0.01;
	R[1] =0.01;
	R[2] =0.01;

	// default magnetometer noise
	R[3] = 0.1;
	R[4] = 0.1;
	R[5] = 0.1;

	magDeclination = 0;
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
	dt2 = 1.0f / (2.0f * freq);
}

/**
 * @brief Sets the gyroscope process noise level for the Extended Kalman Filter.
 *
 * @param[in] NoiseX  Gyroscope process noise variance along the X-axis.
 * @param[in] NoiseY  Gyroscope process noise variance along the Y-axis.
 * @param[in] NoiseZ  Gyroscope process noise variance along the Z-axis.
 *
 * @ return None
 */
void ExtendedKalmanFilter::SetGyroNoise(float NoiseX, float NoiseY, float NoiseZ)
{
	SigmaOmega[0] = NoiseX;
	SigmaOmega[1] = NoiseY;
	SigmaOmega[2] = NoiseZ;
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
  * @brief  Sets the diagonal elements of the estimated state covariance matrix.
  * @param[in]  P00  element P[0][0] (variance of first state component).
  * @param[in]  P11  element P[1][1] (variance of second state component).
  * @param[in]  P22  element P[2][2] (variance of third state component).
  * @param[in]  P33  element P[3][3] (variance of fourth state component).
  *
  * @return None
  *
  * @note This function updates only the diagonal elements of matrix P.
  *       Off-diagonal elements remain unchanged.
  */
void ExtendedKalmanFilter::SetP(float P00, float P11, float P22, float P33)
{
	P[0][0] = P00;
	P[1][1] = P11;
	P[2][2] = P22;
	P[3][3] = P33;
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
 */
bool ExtendedKalmanFilter::Run(float ax, float ay, float az, float gx, float gy, float gz, float mx, float my, float mz)
{
	// Normalize accelerometer and magnetometer readings
	float normA = sqrtf(ax * ax + ay * ay + az * az);
	float normM = sqrtf(mx * mx + my * my + mz * mz);

	if (!normM && !normA) return false;

	if (fabs(normM) >= 0.00001f)
	{
		mx = mx / normM;	my = my / normM; mz = mz / normM;
	}
	else
	{
		Run(ax, ay, az, gx, gy, gz);
	}

	if (fabs(normA) >= 0.00001f)
	{
		ax = ax / normA;	ay = ay / normA;	az = az / normA;
	}


    float mzp = mz + 1.0f;
    float mzm = mz - 1.0f;
    float azp = az + 1.0f;
    float azm = az - 1.0f;

    qam.s = 0.25f * ((ax * mx + ay * my + azp * mzp) * q.s + (-ay * mzm + my * azp) * q.x + (ax * mzm - mx * azp) * q.y + (-ax * my + ay * mx) * q.z);
    qam.x = 0.25f * ((ay * mzp - my * azm) * q.s + (ax * mx + ay * my + azm * mzm) * q.x + (ax * my - ay * mx) * q.y + (ax * mzp - mx * azm) * q.z);
    qam.y = 0.25f * ((-ax * mzp + mx * azm) * q.s + (-ax * my + ay * mx) * q.x + (ax * mx + ay * my + azm * mzm) * q.y + (ay * mzp - my * azm) * q.z);
    qam.z = 0.25f * ((ax * my - ay * mx) * q.s + (-ax * mzm + mx * azp) * q.x + (-ay * mzm + my * azp) * q.y + (ax * mx + ay * my + azp * mzp) * q.z);

    qam.Normalise();
	// convert from degrees/s to rads/s
	gx *= 0.01745f;
	gy *= 0.01745f;
	gz *= 0.01745f;

    // ============================================================
    // PREDICTION STEP - Using gyroscope
    // ============================================================

	// Fill F
    F[0][0] = 1.0;       F[0][1] = -dt2 * gx;  	F[0][2] = -dt2 * gy;  	F[0][3] = -dt2 * gz;
    F[1][0] = dt2 * gx;  F[1][1] = 1.0;         F[1][2] = dt2 * gz;   	F[1][3] = F[0][2];
    F[2][0] = dt2 * gy;  F[2][1] = F[0][3]; 	F[2][2] = 1.0;          F[2][3] = F[1][0];
    F[3][0] = F[1][2];  F[3][1] = F[2][0];  	F[3][2] = F[0][1];  	F[3][3] = 1.0;

	qcap.s = q.s + (F[0][1]*q.x + F[0][2]*q.y + F[0][3]*q.z);
	qcap.x = q.x + (F[1][0]*q.s + F[1][2]*q.y + F[0][2]*q.z);
	qcap.y = q.y + (F[2][0]*q.s + F[1][0]*q.z + F[0][3]*q.x);
	qcap.z = q.z + (F[1][2]*q.s + F[2][0]*q.x + F[0][1]*q.y);

	qcap.Normalise();

	// Process Noise Covariance: Pcap = FPFt + Q
	// FP = F * P
    float FP[4][4] = {0};
	// Pcap = F * P
	for (uint8_t i = 0; i < 4; i++)
		for (uint8_t j = 0; j < 4; j++)
			for (uint8_t k = 0; k < 4; k++)
			{
				if (isnan(P[k][j]))
					return false;
				FP[i][j] += F[i][k] * P[k][j];
			}

	// Finally, Pcap = F*P*Ft + Q
	for (uint8_t i = 0; i < 4; i++)
		for (uint8_t j = 0; j < 4; j++)
		{
			Pcap[i][j] = 0;
			for (uint8_t k = 0; k < 4; k++)
				Pcap[i][j] += FP[i][k] * F[j][k];
			if (isnan(P[i][j]))
				return false;
		}

	Q[0][0] = q.x*SigmaOmega[0]*q.x+q.y*SigmaOmega[1]*q.y+q.z*SigmaOmega[2]*q.z;
	Q[0][1] = q.x*SigmaOmega[0]*-q.s+q.y*SigmaOmega[1]*-q.z+q.z*SigmaOmega[2]*q.y;
	Q[0][2] = q.x*SigmaOmega[0]*q.z+q.y*SigmaOmega[1]*-q.s+q.z*SigmaOmega[2]*-q.x;
	Q[0][3] = q.x*SigmaOmega[0]*-q.y+q.y*SigmaOmega[1]*q.x+q.z*SigmaOmega[2]*-q.s;
	Q[1][0] = -q.s*SigmaOmega[0]*q.x-q.z*SigmaOmega[1]*q.y+q.y*SigmaOmega[2]*q.z;
	Q[1][1] = -q.s*SigmaOmega[0]*-q.s-q.z*SigmaOmega[1]*-q.z+q.y*SigmaOmega[2]*q.y;
	Q[1][2] = -q.s*SigmaOmega[0]*q.z-q.z*SigmaOmega[1]*-q.s+q.y*SigmaOmega[2]*-q.x;
	Q[1][3] = -q.s*SigmaOmega[0]*-q.y-q.z*SigmaOmega[1]*q.x+q.y*SigmaOmega[2]*-q.s;
	Q[2][0] = q.z*SigmaOmega[0]*q.x-q.s*SigmaOmega[1]*q.y-q.x*SigmaOmega[2]*q.z;
	Q[2][1] = q.z*SigmaOmega[0]*-q.s-q.s*SigmaOmega[1]*-q.z-q.x*SigmaOmega[2]*q.y;
	Q[2][2] = q.z*SigmaOmega[0]*q.z-q.s*SigmaOmega[1]*-q.s-q.x*SigmaOmega[2]*-q.x;
	Q[2][3] = q.z*SigmaOmega[0]*-q.y-q.s*SigmaOmega[1]*q.x-q.x*SigmaOmega[2]*-q.s;
	Q[3][0] = -q.y*SigmaOmega[0]*q.x+q.x*SigmaOmega[1]*q.y-q.s*SigmaOmega[2]*q.z;
	Q[3][1] = -q.y*SigmaOmega[0]*-q.s+q.x*SigmaOmega[1]*-q.z-q.s*SigmaOmega[2]*q.y;
	Q[3][2] = -q.y*SigmaOmega[0]*q.z+q.x*SigmaOmega[1]*-q.s-q.s*SigmaOmega[2]*-q.x;
	Q[3][3] = -q.y*SigmaOmega[0]*-q.y+q.x*SigmaOmega[1]*q.x-q.s*SigmaOmega[2]*-q.s;


	Pcap[0][0] += Q[0][0];	Pcap[0][1] += Q[0][1]; Pcap[0][2] += Q[0][2]; Pcap[0][3] += Q[0][3];
	Pcap[1][0] += Q[1][0];	Pcap[1][1] += Q[1][1]; Pcap[1][2] += Q[1][2]; Pcap[1][3] += Q[1][3];
	Pcap[2][0] += Q[2][0];	Pcap[2][1] += Q[2][1]; Pcap[2][2] += Q[2][2]; Pcap[2][3] += Q[2][3];
	Pcap[3][0] += Q[3][0];	Pcap[3][1] += Q[3][1]; Pcap[3][2] += Q[3][2]; Pcap[3][3] += Q[3][3];

    // ============================================================
    // CORRECTION MODEL
    // ============================================================

	// Compute J

    J[0][0] = 0.25f * (mx*q.s - my*q.z + q.y*mzm);
    J[0][1] = 0.25f * (mx*q.z + my*q.s - q.x*mzm);
    J[0][2] = 0.25f * (-mx*q.y + my*q.x + q.s*mzp);
    J[0][3] = 0.25f * (ax*q.s + ay*q.z - q.y*azp);
    J[0][4] = 0.25f * (-ax*q.z + ay*q.s + q.x*azp);
    J[0][5] = 0.25f * (ax*q.y - ay*q.x + q.s*azp);

    J[1][0] = 0.25f * (mx*q.x + my*q.y + q.z*mzp);
    J[1][1] = 0.25f * (-mx*q.y + my*q.x + q.s*mzp);
    J[1][2] = 0.25f * (-mx*q.z - my*q.s + q.x*mzm);
    J[1][3] = 0.25f * (ax*q.x - ay*q.y - q.z*azm);
    J[1][4] = 0.25f * (ax*q.y + ay*q.x - q.s*azm);
    J[1][5] = 0.25f * (ax*q.z + ay*q.s + q.x*azm);

    J[2][0] = 0.25f * (mx*q.y - my*q.x - q.s*mzp);
    J[2][1] = 0.25f * (mx*q.x + my*q.y + q.z*mzp);
    J[2][2] = 0.25f * (mx*q.s - my*q.z + q.y*mzm);
    J[2][3] = 0.25f * (ax*q.y + ay*q.x + q.s*azm);
    J[2][4] = 0.25f * (-ax*q.x + ay*q.y - q.z*azm);
    J[2][5] = 0.25f * (-ax*q.s + ay*q.z + q.y*azm);

    J[3][0] = 0.25f * (mx*q.z + my*q.s - q.x*mzm);
    J[3][1] = 0.25f * (-mx*q.s + my*q.z - q.y*mzm);
    J[3][2] = 0.25f * (mx*q.x + my*q.y + q.z*mzp);
    J[3][3] = 0.25f * (ax*q.z - ay*q.s + q.x*azp);
    J[3][4] = 0.25f * (ax*q.s + ay*q.z + q.y*azp);
    J[3][5] = 0.25f * (-ax*q.x - ay*q.y + q.z*azp);

	JRJt[0][0] = J[0][0]*R[0]*J[0][0]+J[0][1]*R[1]*J[0][1]+J[0][2]*R[2]*J[0][2]+J[0][3]*R[3]*J[0][3]+J[0][4]*R[4]*J[0][4]+J[0][5]*R[5]*J[0][5];
	JRJt[0][1] = J[0][0]*R[0]*J[1][0]+J[0][1]*R[1]*J[1][1]+J[0][2]*R[2]*J[1][2]+J[0][3]*R[3]*J[1][3]+J[0][4]*R[4]*J[1][4]+J[0][5]*R[5]*J[1][5];
	JRJt[0][2] = J[0][0]*R[0]*J[2][0]+J[0][1]*R[1]*J[2][1]+J[0][2]*R[2]*J[2][2]+J[0][3]*R[3]*J[2][3]+J[0][4]*R[4]*J[2][4]+J[0][5]*R[5]*J[2][5];
	JRJt[0][3] = J[0][0]*R[0]*J[3][0]+J[0][1]*R[1]*J[3][1]+J[0][2]*R[2]*J[3][2]+J[0][3]*R[3]*J[3][3]+J[0][4]*R[4]*J[3][4]+J[0][5]*R[5]*J[3][5];

	JRJt[1][0] = J[1][0]*R[0]*J[0][0]+J[1][1]*R[1]*J[0][1]+J[1][2]*R[2]*J[0][2]+J[1][3]*R[3]*J[0][3]+J[1][4]*R[4]*J[0][4]+J[1][5]*R[5]*J[0][5];
	JRJt[1][1] = J[1][0]*R[0]*J[1][0]+J[1][1]*R[1]*J[1][1]+J[1][2]*R[2]*J[1][2]+J[1][3]*R[3]*J[1][3]+J[1][4]*R[4]*J[1][4]+J[1][5]*R[5]*J[1][5];
	JRJt[1][2] = J[1][0]*R[0]*J[2][0]+J[1][1]*R[1]*J[2][1]+J[1][2]*R[2]*J[2][2]+J[1][3]*R[3]*J[2][3]+J[1][4]*R[4]*J[2][4]+J[1][5]*R[5]*J[2][5];
	JRJt[1][3] = J[1][0]*R[0]*J[3][0]+J[1][1]*R[1]*J[3][1]+J[1][2]*R[2]*J[3][2]+J[1][3]*R[3]*J[3][3]+J[1][4]*R[4]*J[3][4]+J[1][5]*R[5]*J[3][5];

	JRJt[2][0] = J[2][0]*R[0]*J[0][0]+J[2][1]*R[1]*J[0][1]+J[2][2]*R[2]*J[0][2]+J[2][3]*R[3]*J[0][3]+J[2][4]*R[4]*J[0][4]+J[2][5]*R[5]*J[0][5];
	JRJt[2][1] = J[2][0]*R[0]*J[1][0]+J[2][1]*R[1]*J[1][1]+J[2][2]*R[2]*J[1][2]+J[2][3]*R[3]*J[1][3]+J[2][4]*R[4]*J[1][4]+J[2][5]*R[5]*J[1][5];
	JRJt[2][2] = J[2][0]*R[0]*J[2][0]+J[2][1]*R[1]*J[2][1]+J[2][2]*R[2]*J[2][2]+J[2][3]*R[3]*J[2][3]+J[2][4]*R[4]*J[2][4]+J[2][5]*R[5]*J[2][5];
	JRJt[2][3] = J[2][0]*R[0]*J[3][0]+J[2][1]*R[1]*J[3][1]+J[2][2]*R[2]*J[3][2]+J[2][3]*R[3]*J[3][3]+J[2][4]*R[4]*J[3][4]+J[2][5]*R[5]*J[3][5];

	JRJt[3][0] = J[3][0]*R[0]*J[0][0]+J[3][1]*R[1]*J[0][1]+J[3][2]*R[2]*J[0][2]+J[3][3]*R[3]*J[0][3]+J[3][4]*R[4]*J[0][4]+J[3][5]*R[5]*J[0][5];
	JRJt[3][1] = J[3][0]*R[0]*J[1][0]+J[3][1]*R[1]*J[1][1]+J[3][2]*R[2]*J[1][2]+J[3][3]*R[3]*J[1][3]+J[3][4]*R[4]*J[1][4]+J[3][5]*R[5]*J[1][5];
	JRJt[3][2] = J[3][0]*R[0]*J[2][0]+J[3][1]*R[1]*J[2][1]+J[3][2]*R[2]*J[2][2]+J[3][3]*R[3]*J[2][3]+J[3][4]*R[4]*J[2][4]+J[3][5]*R[5]*J[2][5];
	JRJt[3][3] = J[3][0]*R[0]*J[3][0]+J[3][1]*R[1]*J[3][1]+J[3][2]*R[2]*J[3][2]+J[3][3]*R[3]*J[3][3]+J[3][4]*R[4]*J[3][4]+J[3][5]*R[5]*J[3][5];


    // ============================================================
    // KALMAN GAIN: Kg = Pcap * S^-1 (4x4 matrix)
    // ============================================================

	float S[4][4] = { 0 }, Sinv[4][4] = { 0 };
	// here S = Pcap + JRJt
	for (uint8_t i = 0; i < 4; ++i)
		for (uint8_t j = 0; j < 4; ++j)
			S[i][j] = Pcap[i][j] + JRJt[i][j];

    // ============================================================
    // MATRIX INVERSION - Compute S^-1
    // ============================================================

	// calculating inverse of S matrix, using LU decomposition
	// S * S^-1 = I, using LU decomposition S = LU
	for (uint8_t p = 0; p < 4; ++p)
	{
		if (fabs(S[p][p]) < 1e-6f)
			return false;

		for (uint8_t i = p + 1; i < 4; ++i)
		{
			//updating the l matrix :l[ur][p]
			S[i][p] /= S[p][p];
			for (uint8_t j = p + 1; j < 4; ++j)
				S[i][j] -= S[i][p] * S[p][j];
		}
	}

	// L * U * S^-1 = I, now let U * S^-1 = Y
	// therefore, L * Y = I
	for (uint8_t col = 0; col < 4; ++col)
	{
		float Y[4] = {0};
		float Icol[4] = {0};
		Icol[col] = 1.0f;

		// Forward substitution: L * Y = Icol
		float sum = 0;
		int8_t i = 0, j = 0;
		for (i = 0; i < 4; ++i)
		{
			sum = 0;
			for (j = 0; j < i; ++j)
				sum += S[i][j] * Y[j];
			Y[i] = Icol[i] - sum;
		}

		// Back substitution: U * X = Y, here X = one column of S^-1 (Sinv)
		for (i = 3; i >= 0; --i)
		{
			sum = 0;
			for (j = i + 1; j < 4; ++j)
				sum += S[i][j] * Sinv[j][col];

			if (fabs(S[i][i]) < 1e-6f)
				return false;

			Sinv[i][col] = (Y[i] - sum) / S[i][i];
		}
	}


 	// Kg = Pcap * S^-1(Sinv)
 	// here, Kg = Pcap * Sinv
 	for (uint8_t m = 0; m < 4; m++)
 		for (uint8_t p = 0; p < 4; p++)
 		{
 			Kg[m][p] = 0;
 			for (uint8_t n = 0; n < 4; n++)
 				Kg[m][p] += Pcap[m][n] * Sinv[n][p];
			if (isnan(Kg[m][p]))
				return false;
 		}

    // ============================================================
    // STATE CORRECTION: q = qcap + Kg (qam - qcap)
    // ============================================================

	// Correction: q = qcap + Kg (qam - qcap), here v = qam - q
 	float v[4];
 	v[0] = qam.s - q.s;
 	v[1] = qam.x - q.x;
 	v[2] = qam.y - q.y;
 	v[3] = qam.z - q.z;

 	q.s = qcap.s + (Kg[0][0]*v[0]+Kg[0][1]*v[1]+Kg[0][2]*v[2]+Kg[0][3]*v[3]);
 	q.x = qcap.x + (Kg[1][0]*v[0]+Kg[1][1]*v[1]+Kg[1][2]*v[2]+Kg[1][3]*v[3]);
 	q.y = qcap.y + (Kg[2][0]*v[0]+Kg[2][1]*v[1]+Kg[2][2]*v[2]+Kg[2][3]*v[3]);
	q.z = qcap.z + (Kg[3][0]*v[0]+Kg[3][1]*v[1]+Kg[3][2]*v[2]+Kg[3][3]*v[3]);

    // ============================================================
    // COVARIANCE UPDATE: P = (I4 - Kg) * Pcap
    // ============================================================

 	// P = (I4 - Kg) * Pcap
 	float _1mKg00 = 1.0f - Kg[0][0];
 	float _1mKg11 = 1.0f - Kg[1][1];
 	float _1mKg22 = 1.0f - Kg[2][2];
 	float _1mKg33 = 1.0f - Kg[3][3];

 	// Now, P = I_Kg * Pcap
 	P[0][0] = _1mKg00*Pcap[0][0]-Kg[0][1]*Pcap[1][0]-Kg[0][2]*Pcap[2][0]-Kg[0][3]*Pcap[3][0];
 	P[0][1] = _1mKg00*Pcap[0][1]-Kg[0][1]*Pcap[1][1]-Kg[0][2]*Pcap[2][1]-Kg[0][3]*Pcap[3][1];
 	P[0][2] = _1mKg00*Pcap[0][2]-Kg[0][1]*Pcap[1][2]-Kg[0][2]*Pcap[2][2]-Kg[0][3]*Pcap[3][2];
 	P[0][3] = _1mKg00*Pcap[0][3]-Kg[0][1]*Pcap[1][3]-Kg[0][2]*Pcap[2][3]-Kg[0][3]*Pcap[3][3];
 	P[1][0] = -Kg[1][0]*Pcap[0][0]+_1mKg11*Pcap[1][0]-Kg[1][2]*Pcap[2][0]-Kg[1][3]*Pcap[3][0];
 	P[1][1] = -Kg[1][0]*Pcap[0][1]+_1mKg11*Pcap[1][1]-Kg[1][2]*Pcap[2][1]-Kg[1][3]*Pcap[3][1];
 	P[1][2] = -Kg[1][0]*Pcap[0][2]+_1mKg11*Pcap[1][2]-Kg[1][2]*Pcap[2][2]-Kg[1][3]*Pcap[3][2];
 	P[1][3] = -Kg[1][0]*Pcap[0][3]+_1mKg11*Pcap[1][3]-Kg[1][2]*Pcap[2][3]-Kg[1][3]*Pcap[3][3];
 	P[2][0] = -Kg[2][0]*Pcap[0][0]-Kg[2][1]*Pcap[1][0]+_1mKg22*Pcap[2][0]-Kg[2][3]*Pcap[3][0];
 	P[2][1] = -Kg[2][0]*Pcap[0][1]-Kg[2][1]*Pcap[1][1]+_1mKg22*Pcap[2][1]-Kg[2][3]*Pcap[3][1];
 	P[2][2] = -Kg[2][0]*Pcap[0][2]-Kg[2][1]*Pcap[1][2]+_1mKg22*Pcap[2][2]-Kg[2][3]*Pcap[3][2];
 	P[2][3] = -Kg[2][0]*Pcap[0][3]-Kg[2][1]*Pcap[1][3]+_1mKg22*Pcap[2][3]-Kg[2][3]*Pcap[3][3];
 	P[3][0] = -Kg[3][0]*Pcap[0][0]-Kg[3][1]*Pcap[1][0]-Kg[3][2]*Pcap[2][0]+_1mKg33*Pcap[3][0];
 	P[3][1] = -Kg[3][0]*Pcap[0][1]-Kg[3][1]*Pcap[1][1]-Kg[3][2]*Pcap[2][1]+_1mKg33*Pcap[3][1];
 	P[3][2] = -Kg[3][0]*Pcap[0][2]-Kg[3][1]*Pcap[1][2]-Kg[3][2]*Pcap[2][2]+_1mKg33*Pcap[3][2];
 	P[3][3] = -Kg[3][0]*Pcap[0][3]-Kg[3][1]*Pcap[1][3]-Kg[3][2]*Pcap[2][3]+_1mKg33*Pcap[3][3];

 	q.Normalise();
 	return true;
}

bool ExtendedKalmanFilter::Run(float ax, float ay, float az, float gx, float gy, float gz)
{
    // Normalize accelerometer reading
    float normA = sqrtf(ax * ax + ay * ay + az * az);
    if (fabs(normA) < 0.00001f)
        return false;

    ax /= normA;
    ay /= normA;
    az /= normA;

    // Compute measurement quaternion qam using accelerometer only
    // Here, assume qam aligns body z-axis with gravity vector (ax, ay, az)
    float azp = az + 1.0f;
    float azm = az - 1.0f;

    qam.s = 0.25 * ( azp * q.s  +  ay * q.x  -  ax * q.y );
    qam.x = 0.25 * ( ay * q.s  -  azm * q.x  +  ax * q.z );
    qam.y = 0.25 * ( -ax * q.s  -  azm * q.y  +  ay * q.z );
    qam.z = 0.25 * ( ax * q.x  +  ay * q.y  +  azp * q.z );

    qam.Normalise();

    // Convert gyro from degrees/s to radians/s
    gx *= 0.01745f;
    gy *= 0.01745f;
    gz *= 0.01745f;

    // ============================================================
    // PREDICTION STEP - Using gyroscope
    // ============================================================

   	// Fill F
    F[0][0] = 1.0;       F[0][1] = -dt2 * gx;  	F[0][2] = -dt2 * gy;  	F[0][3] = -dt2 * gz;
    F[1][0] = dt2 * gx;  F[1][1] = 1.0;         F[1][2] = dt2 * gz;   	F[1][3] = F[0][2];
    F[2][0] = dt2 * gy;  F[2][1] = F[0][3]; 	F[2][2] = 1.0;          F[2][3] = F[1][0];
    F[3][0] = F[1][2];  F[3][1] = F[2][0];  	F[3][2] = F[0][1];  	F[3][3] = 1.0;

	qcap.s = q.s + (F[0][1]*q.x + F[0][2]*q.y + F[0][3]*q.z);
	qcap.x = q.x + (F[1][0]*q.s + F[1][2]*q.y + F[0][2]*q.z);
	qcap.y = q.y + (F[2][0]*q.s + F[1][0]*q.z + F[0][3]*q.x);
	qcap.z = q.z + (F[1][2]*q.s + F[2][0]*q.x + F[0][1]*q.y);

   	qcap.Normalise();

   	// Process Noise Covariance: Pcap = FPFt + Q
   	// FP = F * P
       float FP[4][4] = {0};
   	// Pcap = F * P
   	for (uint8_t i = 0; i < 4; i++)
   		for (uint8_t j = 0; j < 4; j++)
   			for (uint8_t k = 0; k < 4; k++)
   			{
   				if (isnan(P[k][j]))
   					return false;
   				FP[i][j] += F[i][k] * P[k][j];
   			}

   	// Finally, Pcap = F*P*Ft + Q
   	for (uint8_t i = 0; i < 4; i++)
   		for (uint8_t j = 0; j < 4; j++)
   		{
   			Pcap[i][j] = 0;
   			for (uint8_t k = 0; k < 4; k++)
   				Pcap[i][j] += FP[i][k] * F[j][k];
   			if (isnan(P[i][j]))
   				return false;
   		}

   	Q[0][0] = q.x*SigmaOmega[0]*q.x+q.y*SigmaOmega[1]*q.y+q.z*SigmaOmega[2]*q.z;
   	Q[0][1] = q.x*SigmaOmega[0]*-q.s+q.y*SigmaOmega[1]*-q.z+q.z*SigmaOmega[2]*q.y;
   	Q[0][2] = q.x*SigmaOmega[0]*q.z+q.y*SigmaOmega[1]*-q.s+q.z*SigmaOmega[2]*-q.x;
   	Q[0][3] = q.x*SigmaOmega[0]*-q.y+q.y*SigmaOmega[1]*q.x+q.z*SigmaOmega[2]*-q.s;
   	Q[1][0] = -q.s*SigmaOmega[0]*q.x-q.z*SigmaOmega[1]*q.y+q.y*SigmaOmega[2]*q.z;
   	Q[1][1] = -q.s*SigmaOmega[0]*-q.s-q.z*SigmaOmega[1]*-q.z+q.y*SigmaOmega[2]*q.y;
   	Q[1][2] = -q.s*SigmaOmega[0]*q.z-q.z*SigmaOmega[1]*-q.s+q.y*SigmaOmega[2]*-q.x;
   	Q[1][3] = -q.s*SigmaOmega[0]*-q.y-q.z*SigmaOmega[1]*q.x+q.y*SigmaOmega[2]*-q.s;
   	Q[2][0] = q.z*SigmaOmega[0]*q.x-q.s*SigmaOmega[1]*q.y-q.x*SigmaOmega[2]*q.z;
   	Q[2][1] = q.z*SigmaOmega[0]*-q.s-q.s*SigmaOmega[1]*-q.z-q.x*SigmaOmega[2]*q.y;
   	Q[2][2] = q.z*SigmaOmega[0]*q.z-q.s*SigmaOmega[1]*-q.s-q.x*SigmaOmega[2]*-q.x;
   	Q[2][3] = q.z*SigmaOmega[0]*-q.y-q.s*SigmaOmega[1]*q.x-q.x*SigmaOmega[2]*-q.s;
   	Q[3][0] = -q.y*SigmaOmega[0]*q.x+q.x*SigmaOmega[1]*q.y-q.s*SigmaOmega[2]*q.z;
   	Q[3][1] = -q.y*SigmaOmega[0]*-q.s+q.x*SigmaOmega[1]*-q.z-q.s*SigmaOmega[2]*q.y;
   	Q[3][2] = -q.y*SigmaOmega[0]*q.z+q.x*SigmaOmega[1]*-q.s-q.s*SigmaOmega[2]*-q.x;
   	Q[3][3] = -q.y*SigmaOmega[0]*-q.y+q.x*SigmaOmega[1]*q.x-q.s*SigmaOmega[2]*-q.s;


   	Pcap[0][0] += Q[0][0];	Pcap[0][1] += Q[0][1]; Pcap[0][2] += Q[0][2]; Pcap[0][3] += Q[0][3];
   	Pcap[1][0] += Q[1][0];	Pcap[1][1] += Q[1][1]; Pcap[1][2] += Q[1][2]; Pcap[1][3] += Q[1][3];
   	Pcap[2][0] += Q[2][0];	Pcap[2][1] += Q[2][1]; Pcap[2][2] += Q[2][2]; Pcap[2][3] += Q[2][3];
   	Pcap[3][0] += Q[3][0];	Pcap[3][1] += Q[3][1]; Pcap[3][2] += Q[3][2]; Pcap[3][3] += Q[3][3];

    // ============================================================
    // CORRECTION MODEL
    // ============================================================

   	// without the magnetometer J becomes a 4x3 matrix
   	J[0][0] = 0.25f *  -q.y;
   	J[0][1] = 0.25f * q.x;
   	J[0][2] = 0.25f * q.s;

   	J[1][0] = 0.25f * q.z;
   	J[1][1] = 0.25f * q.s;
   	J[1][2] = 0.25f * -q.x;

   	J[2][0] = 0.25f * -q.s;
   	J[2][1] = 0.25f * q.z;
   	J[2][2] = 0.25f * -q.y;

   	J[3][0] = 0.25f * q.x;
   	J[3][1] = 0.25f * q.y;
   	J[3][2] = 0.25f * q.z;


	JRJt[0][0] = J[0][0]*R[0]*J[0][0]+J[0][1]*R[1]*J[0][1]+J[0][2]*R[2]*J[0][2];
	JRJt[0][1] = J[0][0]*R[0]*J[1][0]+J[0][1]*R[1]*J[1][1]+J[0][2]*R[2]*J[1][2];
	JRJt[0][2] = J[0][0]*R[0]*J[2][0]+J[0][1]*R[1]*J[2][1]+J[0][2]*R[2]*J[2][2];

	JRJt[1][0] = J[1][0]*R[0]*J[0][0]+J[1][1]*R[1]*J[0][1]+J[1][2]*R[2]*J[0][2];
	JRJt[1][1] = J[1][0]*R[0]*J[1][0]+J[1][1]*R[1]*J[1][1]+J[1][2]*R[2]*J[1][2];
	JRJt[1][2] = J[1][0]*R[0]*J[2][0]+J[1][1]*R[1]*J[2][1]+J[1][2]*R[2]*J[2][2];

	JRJt[2][0] = J[2][0]*R[0]*J[0][0]+J[2][1]*R[1]*J[0][1]+J[2][2]*R[2]*J[0][2];
	JRJt[2][1] = J[2][0]*R[0]*J[1][0]+J[2][1]*R[1]*J[1][1]+J[2][2]*R[2]*J[1][2];
	JRJt[2][2] = J[2][0]*R[0]*J[2][0]+J[2][1]*R[1]*J[2][1]+J[2][2]*R[2]*J[2][2];

	JRJt[3][0] = J[3][0]*R[0]*J[0][0]+J[3][1]*R[1]*J[0][1]+J[3][2]*R[2]*J[0][2];
	JRJt[3][1] = J[3][0]*R[0]*J[1][0]+J[3][1]*R[1]*J[1][1]+J[3][2]*R[2]*J[1][2];
	JRJt[3][2] = J[3][0]*R[0]*J[2][0]+J[3][1]*R[1]*J[2][1]+J[3][2]*R[2]*J[2][2];


	// ============================================================
	// KALMAN GAIN: Kg = Pcap * S^-1 (4x4 matrix)
	// ============================================================

	float S[4][4] = { 0 }, Sinv[4][4] = { 0 };
	// here S = Pcap + JRJt
	for (uint8_t i = 0; i < 4; ++i)
		for (uint8_t j = 0; j < 4; ++j)
			S[i][j] = Pcap[i][j] + JRJt[i][j];

	// ============================================================
	// MATRIX INVERSION - Compute S^-1
	// ============================================================

	// calculating inverse of S matrix, using LU decomposition
	// S * S^-1 = I, using LU decomposition S = LU
	for (uint8_t p = 0; p < 4; ++p)
	{
		if (fabs(S[p][p]) < 1e-6f)
			return false;

		for (uint8_t i = p + 1; i < 4; ++i)
		{
			//updating the l matrix :l[ur][p]
			S[i][p] /= S[p][p];
			for (uint8_t j = p + 1; j < 4; ++j)
				S[i][j] -= S[i][p] * S[p][j];
		}
	}

	// L * U * S^-1 = I, now let U * S^-1 = Y
	// therefore, L * Y = I
	for (uint8_t col = 0; col < 4; ++col)
	{
		float Y[4] = {0};
		float Icol[4] = {0};
		Icol[col] = 1.0f;

		// Forward substitution: L * Y = Icol
		float sum = 0;
		int8_t i = 0, j = 0;
		for (i = 0; i < 4; ++i)
		{
			sum = 0;
			for (j = 0; j < i; ++j)
				sum += S[i][j] * Y[j];
			Y[i] = Icol[i] - sum;
		}

		// Back substitution: U * X = Y, here X = one column of S^-1 (Sinv)
		for (i = 3; i >= 0; --i)
		{
			sum = 0;
			for (j = i + 1; j < 4; ++j)
				sum += S[i][j] * Sinv[j][col];

			if (fabs(S[i][i]) < 1e-6f)
				return false;

			Sinv[i][col] = (Y[i] - sum) / S[i][i];
		}
	}


	// Kg = Pcap * S^-1(Sinv)
	// here, Kg = Pcap * Sinv
	for (uint8_t m = 0; m < 4; m++)
		for (uint8_t p = 0; p < 4; p++)
		{
			Kg[m][p] = 0;
			for (uint8_t n = 0; n < 4; n++)
				Kg[m][p] += Pcap[m][n] * Sinv[n][p];
			if (isnan(Kg[m][p]))
				return false;
		}

	// ============================================================
	// STATE CORRECTION: q = qcap + Kg (qam - qcap)
	// ============================================================

	// Correction: q = qcap + Kg (qam - qcap), here v = qam - q
	float v[4];
	v[0] = qam.s - q.s;
	v[1] = qam.x - q.x;
	v[2] = qam.y - q.y;
	v[3] = qam.z - q.z;

	q.s = qcap.s + (Kg[0][0]*v[0]+Kg[0][1]*v[1]+Kg[0][2]*v[2]+Kg[0][3]*v[3]);
	q.x = qcap.x + (Kg[1][0]*v[0]+Kg[1][1]*v[1]+Kg[1][2]*v[2]+Kg[1][3]*v[3]);
	q.y = qcap.y + (Kg[2][0]*v[0]+Kg[2][1]*v[1]+Kg[2][2]*v[2]+Kg[2][3]*v[3]);
	q.z = qcap.z + (Kg[3][0]*v[0]+Kg[3][1]*v[1]+Kg[3][2]*v[2]+Kg[3][3]*v[3]);

	// ============================================================
	// COVARIANCE UPDATE: P = (I4 - Kg) * Pcap
	// ============================================================

	// P = (I4 - Kg) * Pcap
	float _1mKg00 = 1.0f - Kg[0][0];
	float _1mKg11 = 1.0f - Kg[1][1];
	float _1mKg22 = 1.0f - Kg[2][2];
	float _1mKg33 = 1.0f - Kg[3][3];

	// Now, P = I_Kg * Pcap
	P[0][0] = _1mKg00*Pcap[0][0]-Kg[0][1]*Pcap[1][0]-Kg[0][2]*Pcap[2][0]-Kg[0][3]*Pcap[3][0];
	P[0][1] = _1mKg00*Pcap[0][1]-Kg[0][1]*Pcap[1][1]-Kg[0][2]*Pcap[2][1]-Kg[0][3]*Pcap[3][1];
	P[0][2] = _1mKg00*Pcap[0][2]-Kg[0][1]*Pcap[1][2]-Kg[0][2]*Pcap[2][2]-Kg[0][3]*Pcap[3][2];
	P[0][3] = _1mKg00*Pcap[0][3]-Kg[0][1]*Pcap[1][3]-Kg[0][2]*Pcap[2][3]-Kg[0][3]*Pcap[3][3];
	P[1][0] = -Kg[1][0]*Pcap[0][0]+_1mKg11*Pcap[1][0]-Kg[1][2]*Pcap[2][0]-Kg[1][3]*Pcap[3][0];
	P[1][1] = -Kg[1][0]*Pcap[0][1]+_1mKg11*Pcap[1][1]-Kg[1][2]*Pcap[2][1]-Kg[1][3]*Pcap[3][1];
	P[1][2] = -Kg[1][0]*Pcap[0][2]+_1mKg11*Pcap[1][2]-Kg[1][2]*Pcap[2][2]-Kg[1][3]*Pcap[3][2];
	P[1][3] = -Kg[1][0]*Pcap[0][3]+_1mKg11*Pcap[1][3]-Kg[1][2]*Pcap[2][3]-Kg[1][3]*Pcap[3][3];
	P[2][0] = -Kg[2][0]*Pcap[0][0]-Kg[2][1]*Pcap[1][0]+_1mKg22*Pcap[2][0]-Kg[2][3]*Pcap[3][0];
	P[2][1] = -Kg[2][0]*Pcap[0][1]-Kg[2][1]*Pcap[1][1]+_1mKg22*Pcap[2][1]-Kg[2][3]*Pcap[3][1];
	P[2][2] = -Kg[2][0]*Pcap[0][2]-Kg[2][1]*Pcap[1][2]+_1mKg22*Pcap[2][2]-Kg[2][3]*Pcap[3][2];
	P[2][3] = -Kg[2][0]*Pcap[0][3]-Kg[2][1]*Pcap[1][3]+_1mKg22*Pcap[2][3]-Kg[2][3]*Pcap[3][3];
	P[3][0] = -Kg[3][0]*Pcap[0][0]-Kg[3][1]*Pcap[1][0]-Kg[3][2]*Pcap[2][0]+_1mKg33*Pcap[3][0];
	P[3][1] = -Kg[3][0]*Pcap[0][1]-Kg[3][1]*Pcap[1][1]-Kg[3][2]*Pcap[2][1]+_1mKg33*Pcap[3][1];
	P[3][2] = -Kg[3][0]*Pcap[0][2]-Kg[3][1]*Pcap[1][2]-Kg[3][2]*Pcap[2][2]+_1mKg33*Pcap[3][2];
	P[3][3] = -Kg[3][0]*Pcap[0][3]-Kg[3][1]*Pcap[1][3]-Kg[3][2]*Pcap[2][3]+_1mKg33*Pcap[3][3];

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
void ExtendedKalmanFilter::GetOrientation(Quaternion& qState)
{
	qState.s = q.s;
	qState.x = q.x;
	qState.y = q.y;
	qState.z = q.z;
}
