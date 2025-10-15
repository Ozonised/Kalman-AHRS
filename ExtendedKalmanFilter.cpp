#include "ExtendedKalmanFilter.hpp"

ExtendedKalmanFilter::ExtendedKalmanFilter()
{
	q.s = 1;
	qcap.s = 1;

	memset(P, 0, sizeof(P));
	P[0][0] = 0.03;
	P[1][1] = 0.03;
	P[2][2] = 0.03;
	P[3][3] = 0.03;

	memset(a, 0, sizeof(a));
	memset(m, 0, sizeof(m));
	memset(U, 0, sizeof(U));
	// default gyroscoppe noise
	SigmaOmega = 0.01;

	// default accelerometer noise
	R[0] =0.05;
	R[1] =0.05;
	R[2] =0.05;

	// default magnetometer noise
	R[3] = 0.8;
	R[4] = 0.8;
	R[5] = 0.8;

	magDeclination = 0;
	rx = 1;
	rz = 0;
}

void ExtendedKalmanFilter::SetSampleTime(float t)
{
	dt = 1.0f / t;
	dt2 = dt / 2;
}

void ExtendedKalmanFilter::SetMagneticDip(float degrees)
{
	float rad = degrees *  0.01745f;
	rx = cos(rad);
	rz = sin(rad);
}

void ExtendedKalmanFilter::SetGyroNoise(float Noise)
{
	SigmaOmega = Noise;
}

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

bool ExtendedKalmanFilter::Run(float ax, float ay, float az, float gx, float gy, float gz, float mx, float my, float mz)
{
	// Normalize accelerometer and gyroscope readings
	float normA = sqrtf(ax * ax + ay * ay + az * az);
	float normM = sqrtf(mx * mx + my * my + mz * mz);

	if (normA)
	{
		a[0] = ax / normA;	a[1] = ay / normA;	a[2] = az / normA;
	}

	if (normM)
	{
		m[0] = mx / normM;	m[1] = my / normM; m[2] = mz / normM;
	}

	U[0] = gx; U[1] = gy; U[2] = gz;


    // ============================================================
    // PREDICTION STEP - Using gyroscope
    // ============================================================

	qcap.s = q.s - dt2 * (U[0]*q.x + U[1]*q.y + U[2]*q.z);
	qcap.x = q.x + dt2 * (U[0]*q.s - U[1]*q.z + U[2]*q.y);
	qcap.y = q.y + dt2 * (U[0]*q.z + U[1]*q.s - U[2]*q.x);
	qcap.z = q.z + dt2 * (-U[0]*q.y + U[1]*q.x + U[2]*q.s);

	qcap.Normalise();

	// Fill F
    F[0][0] = 1.0;       F[0][1] = -dt2 * U[0];  F[0][2] = -dt2 * U[1];  F[0][3] = -dt2 * U[2];
    F[1][0] =  dt2 * U[0]; F[1][1] = 1.0;        F[1][2] =  dt2 * U[2];  F[1][3] = -dt2 * U[1];
    F[2][0] =  dt2 * U[1]; F[2][1] = -dt2 * U[2]; F[2][2] = 1.0;         F[2][3] =  dt2 * U[0];
    F[3][0] =  dt2 * U[2]; F[3][1] =  dt2 * U[1]; F[3][2] = -dt2 * U[0];  F[3][3] = 1.0;

	// Process Noise Covariance: Pcap = FPFt + Q
	// FP = F * P
    float FP[4][4] = {0}, FPFt[4][4] = {0};


	// Pcap = F * P

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			for (int k = 0; k < 4; k++)
			{
				if (isnan(P[k][j]))
					return false;
				FP[i][j] += F[i][k] * P[k][j];
			}

	// Finally, Pcap = F*P*Ft + Q
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
		{
			FPFt[i][j] = 0;
			for (int k = 0; k < 4; k++)
				FPFt[i][j] += FP[i][k] * F[j][k];
			Pcap[i][j] = FPFt[i][j];
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
	v[0] = a[0] - acap[0];
	v[1] = a[1] - acap[1];
	v[2] = a[2] - acap[2];

	if (normM)
	{
		v[3] = m[0] - mcap[0];
		v[4] = m[1] - mcap[1];
		v[5] = m[2] - mcap[2];
	} else
	{
        // If no magnetometer, zero out these residuals to ignore magnetometer update
        v[3] = 0.0f;
        v[4] = 0.0f;
        v[5] = 0.0f;
	}

	// S = H(qcap)*Pcap*H(qcap)^T + R
	// H(qcap) = Jacobian of h(qcap)
	Hqcap[0][0] = 2 * qcap.y;	Hqcap[0][1] = 2 * qcap.z;	Hqcap[0][2] = 2 * qcap.s;	Hqcap[0][3] = 2 * qcap.x;
	Hqcap[1][0] = -2 * qcap.x;	Hqcap[1][1] = -2 * qcap.s;	Hqcap[1][2] = 2 * qcap.z;	Hqcap[1][3] = 2 * qcap.y;
	Hqcap[2][0] = 2 * qcap.s;	Hqcap[2][1] = 2 * -qcap.x;	Hqcap[2][2] = 2 * -qcap.y;	Hqcap[2][3] = 2 * qcap.z;

	if (normM)
	{
		Hqcap[3][0] = 2.0f * (rx*qcap.s - rz*qcap.y);  Hqcap[3][1] = 2.0f * (rx*qcap.x + rz*qcap.z);  Hqcap[3][2] = 2.0f * (-rx*qcap.y - rz*qcap.s); Hqcap[3][3] = 2.0f * (rx*qcap.z - rz*qcap.x);
		Hqcap[4][0] = 2.0f * (rx*qcap.z + rz*qcap.x);  Hqcap[4][1] = 2.0f * (rx*qcap.y + rz*qcap.s);  Hqcap[4][2] = 2.0f * (rx*qcap.x - rz*qcap.y);  Hqcap[4][3] = 2.0f * (rx*qcap.s + rz*qcap.z);
		Hqcap[5][0] = 2.0f * (-rx*qcap.y + rz*qcap.s); Hqcap[5][1] = 2.0f * (rx*qcap.z - rz*qcap.x);  Hqcap[5][2] = 2.0f * (rx*qcap.s + rz*qcap.y);  Hqcap[5][3] = 2.0f * (rx*qcap.y - rz*qcap.s);
	} else
	{
        // Zero out magnetometer Jacobian rows to ignore magnetometer update
        for (int i = 3; i < 6; i++)
            for (int j = 0; j < 4; j++)
                Hqcap[i][j] = 0.0f;
	}

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
 			sinv[i][col] = Y[i];
 			for (int8_t j = i + 1; j < 6; j++)
 				sinv[i][col] -= s[i][j] * sinv[j][col];

 			// singular pivot
 			if (fabs(s[i][i]) <= 0.00001f)
 				return false;

 			sinv[i][col] /= s[i][i];
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
 				Kg[m][p] += PcapHT[m][n] * sinv[n][p];
			if (isnan(Kg[m][p]))
				return false;
 		}

    // ============================================================
    // STATE CORRECTION: q = qcap + Kg * v
    // ============================================================

	// Correction: q = qcap + Kg(z - h(q))
 	// q = qcap + Kg * V
 	// q = Kg * v
 	q.s = Kg[0][0]*v[0]+Kg[0][1]*v[1]+Kg[0][2]*v[2]+Kg[0][3]*v[3]+Kg[0][4]*v[4]+Kg[0][5]*v[5];
 	q.x = Kg[1][0]*v[0]+Kg[1][1]*v[1]+Kg[1][2]*v[2]+Kg[1][3]*v[3]+Kg[1][4]*v[4]+Kg[1][5]*v[5];
 	q.y = Kg[2][0]*v[0]+Kg[2][1]*v[1]+Kg[2][2]*v[2]+Kg[2][3]*v[3]+Kg[2][4]*v[4]+Kg[2][5]*v[5];
 	q.z = Kg[3][0]*v[0]+Kg[3][1]*v[1]+Kg[3][2]*v[2]+Kg[3][3]*v[3]+Kg[3][4]*v[4]+Kg[3][5]*v[5];

 	// q = q + qcap
	q.s += qcap.s;
 	q.x += qcap.x;
 	q.y += qcap.y;
 	q.z	+= qcap.z;

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

void ExtendedKalmanFilter::GetOrientation(Quanternion& qState)
{
	qState.s = q.s;
	qState.x = q.x;
	qState.y = q.y;
	qState.z = q.z;
}
