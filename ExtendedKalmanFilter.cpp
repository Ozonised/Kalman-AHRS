#include "ExtendedKalmanFilter.hpp"

ExtendedKalmanFilter::ExtendedKalmanFilter()
{
	memset(P, 0, sizeof(P));
	P[0][0] = 1;
	P[1][1] = 1;
	P[2][2] = 1;
	P[3][3] = 1;

	// default gyroscoppe noise
	SigmaOmega[0] = 0.3 * 0.3;
	SigmaOmega[1] = 0.3 * 0.3;
	SigmaOmega[2] = 0.3 * 0.3;

	// default accelerometer noise
	R[0] = 0.5 * 0.5;
	R[1] = 0.5 * 0.5;
	R[2] = 0.5 * 0.5;

	// default magnetometer noise
	R[3] = 0.8 * 0.8;
	R[4] = 0.8 * 0.8;
	R[5] = 0.8 * 0.8;
}

void ExtendedKalmanFilter::SetSampleTime(float t)
{
	dt2 = t / 2;
}

void ExtendedKalmanFilter::SetMagneticInclination(float degrees)
{
	float rad = degrees *  0.01745f;
	rx = cos(rad);
	rz = sin(rad);
}

void ExtendedKalmanFilter::SetGyroNoise(float x, float y, float z)
{
	SigmaOmega[0] = x;
	SigmaOmega[1] = y;
	SigmaOmega[2] = z;
}

bool ExtendedKalmanFilter::Run(float ax, float ay, float az, float mz, float gx, float gy, float gz, float mx, float my)
{
	float qcapS_Squared = qcap.s * qcap.s;
	float qcapX_Squared = qcap.x * qcap.x;
	float qcapY_Squared = qcap.y * qcap.y;
	float qcapZ_Squared = qcap.z * qcap.z;

	// Normalize accelerometer and gyroscope readings
	float normA = sqrtf(ax * ax + ay * ay + az * az);
	float normM = sqrtf(mx * mx + my * my + mz * mz);

	a[0] = ax / normA;	a[1] = ay / normA;	a[2] = az / normA;
	m[0] = mx / normM;	m[0] = my / normM; m[0] = mz / normM;

	qcap.s = q.s - dt2 * U[0] * q.x - dt2 * U[1] * q.y - dt2 * U[2] * q.z;
	qcap.x = q.x + dt2 * U[0] * q.s - dt2 * U[1] * q.z + dt2 * U[2] * q.y;
	qcap.y = q.y + dt2 * U[0] * q.z + dt2 * U[1] * q.s - dt2 * U[2] * q.x;
	qcap.z = q.z - dt2 * U[0] * q.y + dt2 * U[1] * q.x + dt2 * U[2] * q.s;

	// Fill F
    F[0][0] = 1.0;       F[0][1] = -dt2 * U[0];  F[0][2] = -dt2 * U[1];  F[0][3] = -dt2 * U[2];
    F[1][0] =  dt2 * U[0]; F[1][1] = 1.0;        F[1][2] =  dt2 * U[2];  F[1][3] = -dt2 * U[1];
    F[2][0] =  dt2 * U[1]; F[2][1] = -dt2 * U[2]; F[2][2] = 1.0;         F[2][3] =  dt2 * U[0];
    F[3][0] =  dt2 * U[2]; F[3][1] =  dt2 * U[1]; F[3][2] = -dt2 * U[0];  F[3][3] = 1.0;

	// Process Noise Covariance: Pcap = FPFt + Q
	// FP = F * P
	FP[0][0] = (F[0][0]*P[0][0])+(F[0][1]*P[1][0])+(F[0][2]*P[2][0])+(F[0][3]*P[3][0]);
	FP[0][1] = (F[0][0]*P[0][1])+(F[0][1]*P[1][1])+(F[0][2]*P[2][1])+(F[0][3]*P[3][1]);
	FP[0][2] = (F[0][0]*P[0][2])+(F[0][1]*P[1][2])+(F[0][2]*P[2][2])+(F[0][3]*P[3][2]);
	FP[0][3] = (F[0][0]*P[0][3])+(F[0][1]*P[1][3])+(F[0][2]*P[2][3])+(F[0][3]*P[3][3]);

	FP[1][0] = (F[1][0]*P[0][0])+(F[1][1]*P[1][0])+(F[1][2]*P[2][0])+(F[1][3]*P[3][0]);
	FP[1][1] = (F[1][0]*P[0][1])+(F[1][1]*P[1][1])+(F[1][2]*P[2][1])+(F[1][3]*P[3][1]);
	FP[1][2] = (F[1][0]*P[0][2])+(F[1][1]*P[1][2])+(F[1][2]*P[2][2])+(F[1][3]*P[3][2]);
	FP[1][3] = (F[1][0]*P[0][3])+(F[1][1]*P[1][3])+(F[1][2]*P[2][3])+(F[1][3]*P[3][3]);

	FP[2][0] = (F[2][0]*P[0][0])+(F[2][1]*P[1][0])+(F[2][2]*P[2][0])+(F[2][3]*P[3][0]);
	FP[2][1] = (F[2][0]*P[0][1])+(F[2][1]*P[1][1])+(F[2][2]*P[2][1])+(F[2][3]*P[3][1]);
	FP[2][2] = (F[2][0]*P[0][2])+(F[2][1]*P[1][2])+(F[2][2]*P[2][2])+(F[2][3]*P[3][2]);
	FP[2][3] = (F[2][0]*P[0][3])+(F[2][1]*P[1][3])+(F[2][2]*P[2][3])+(F[2][3]*P[3][3]);

	FP[3][0] = (F[3][0]*P[0][0])+(F[3][1]*P[1][0])+(F[3][2]*P[2][0])+(F[3][3]*P[3][0]);
	FP[3][1] = (F[3][0]*P[0][1])+(F[3][1]*P[1][1])+(F[3][2]*P[2][1])+(F[3][3]*P[3][1]);
	FP[3][2] = (F[3][0]*P[0][2])+(F[3][1]*P[1][2])+(F[3][2]*P[2][2])+(F[3][3]*P[3][2]);
	FP[3][3] = (F[3][0]*P[0][3])+(F[3][1]*P[1][3])+(F[3][2]*P[2][3])+(F[3][3]*P[3][3]);

//	// FPFt = F*P*Ft
//	FPFt[0][0] = (FP[0][0]*F[0][0])+(FP[0][1]*F[0][1])+(FP[0][2]*F[0][2])+(FP[0][3]*F[0][3]);
//	FPFt[0][1] = (FP[0][0]*F[1][0])+(FP[0][1]*F[1][1])+(FP[0][2]*F[1][2])+(FP[0][3]*F[1][3]);
//	FPFt[0][2] = (FP[0][0]*F[2][0])+(FP[0][1]*F[2][1])+(FP[0][2]*F[2][2])+(FP[0][3]*F[2][3]);
//	FPFt[0][3] = (FP[0][0]*F[3][0])+(FP[0][1]*F[3][1])+(FP[0][2]*F[3][2])+(FP[0][3]*F[3][3]);
//
//	FPFt[1][0] = (FP[1][0]*F[0][0])+(FP[1][1]*F[0][1])+(FP[1][2]*F[0][2])+(FP[1][3]*F[0][3]);
//	FPFt[1][1] = (FP[1][0]*F[1][0])+(FP[1][1]*F[1][1])+(FP[1][2]*F[1][2])+(FP[1][3]*F[1][3]);
//	FPFt[1][2] = (FP[1][0]*F[2][0])+(FP[1][1]*F[2][1])+(FP[1][2]*F[2][2])+(FP[1][3]*F[2][3]);
//	FPFt[1][3] = (FP[1][0]*F[3][0])+(FP[1][1]*F[3][1])+(FP[1][2]*F[3][2])+(FP[1][3]*F[3][3]);
//
//	FPFt[2][0] = (FP[2][0]*F[0][0])+(FP[2][1]*F[0][1])+(FP[2][2]*F[0][2])+(FP[2][3]*F[0][3]);
//	FPFt[2][1] = (FP[2][0]*F[1][0])+(FP[2][1]*F[1][1])+(FP[2][2]*F[1][2])+(FP[2][3]*F[1][3]);
//	FPFt[2][2] = (FP[2][0]*F[2][0])+(FP[2][1]*F[2][1])+(FP[2][2]*F[2][2])+(FP[2][3]*F[2][3]);
//	FPFt[2][3] = (FP[2][0]*F[3][0])+(FP[2][1]*F[3][1])+(FP[2][2]*F[3][2])+(FP[2][3]*F[3][3]);
//
//	FPFt[3][0] = (FP[3][0]*F[0][0])+(FP[3][1]*F[0][1])+(FP[3][2]*F[0][2])+(FP[3][3]*F[0][3]);
//	FPFt[3][1] = (FP[3][0]*F[1][0])+(FP[3][1]*F[1][1])+(FP[3][2]*F[1][2])+(FP[3][3]*F[1][3]);
//	FPFt[3][2] = (FP[3][0]*F[2][0])+(FP[3][1]*F[2][1])+(FP[3][2]*F[2][2])+(FP[3][3]*F[2][3]);
//	FPFt[3][3] = (FP[3][0]*F[3][0])+(FP[3][1]*F[3][1])+(FP[3][2]*F[3][2])+(FP[3][3]*F[3][3]);

	// Q = W * SigmaOmega * Wt
	Q[0][0] = (-q.x*SigmaOmega[0]*-q.x)+(-q.y*SigmaOmega[1]*-q.y)+(-q.z*SigmaOmega[2]*-q.z);
	Q[0][1] = (-q.x*SigmaOmega[0]*q.s)+(-q.y*SigmaOmega[1]*-q.z)+(-q.z*SigmaOmega[2]*q.y);
	Q[0][2] = (-q.x*SigmaOmega[0]*q.z)+(-q.y*SigmaOmega[1]*q.s)+(-q.z*SigmaOmega[2]*-q.x);
	Q[0][3] = (-q.x*SigmaOmega[0]*-q.y)+(-q.y*SigmaOmega[1]*q.s)+(-q.z*SigmaOmega[2]*q.s);

	Q[1][0] = (q.s*SigmaOmega[0]*-q.x)+(-q.z*SigmaOmega[1]*-q.y)+(q.y*SigmaOmega[2]*-q.z);
	Q[1][1] = (q.s*SigmaOmega[0]*q.s)+(-q.z*SigmaOmega[1]*-q.z)+(q.y*SigmaOmega[2]*q.y);
	Q[1][2] = (q.s*SigmaOmega[0]*q.z)+(-q.z*SigmaOmega[1]*q.s)+(q.y*SigmaOmega[2]*-q.x);
	Q[1][3] = (q.s*SigmaOmega[0]*-q.y)+(-q.z*SigmaOmega[1]*q.s)+(q.y*SigmaOmega[2]*q.s);

	Q[2][0] = (q.z*SigmaOmega[0]*-q.x)+(q.s*SigmaOmega[1]*-q.y)+(-q.x*SigmaOmega[2]*-q.z);
	Q[2][1] = (q.z*SigmaOmega[0]*q.s)+(q.s*SigmaOmega[1]*-q.z)+(-q.x*SigmaOmega[2]*q.y);
	Q[2][2] = (q.z*SigmaOmega[0]*q.z)+(q.s*SigmaOmega[1]*q.s)+(-q.x*SigmaOmega[2]*-q.x);
	Q[2][3] = (q.z*SigmaOmega[0]*-q.y)+(q.s*SigmaOmega[1]*q.s)+(-q.x*SigmaOmega[2]*q.s);

	Q[3][0] = (-q.y*SigmaOmega[0]*-q.x)+(q.x*SigmaOmega[1]*-q.y)+(q.s*SigmaOmega[2]*-q.z);
	Q[3][1] = (-q.y*SigmaOmega[0]*q.s)+(q.x*SigmaOmega[1]*-q.z)+(q.s*SigmaOmega[2]*q.y);
	Q[3][2] = (-q.y*SigmaOmega[0]*q.z)+(q.x*SigmaOmega[1]*q.s)+(q.s*SigmaOmega[2]*-q.x);
	Q[3][3] = (-q.y*SigmaOmega[0]*-q.y)+(q.x*SigmaOmega[1]*q.s)+(q.s*SigmaOmega[2]*q.s);

	// Finally, Pcap = F*P*Ft + Q
	Pcap[0][0] = ((FP[0][0]*F[0][0])+(FP[0][1]*F[0][1])+(FP[0][2]*F[0][2])+(FP[0][3]*F[0][3])) + Q[0][0];
	Pcap[0][1] = ((FP[0][0]*F[1][0])+(FP[0][1]*F[1][1])+(FP[0][2]*F[1][2])+(FP[0][3]*F[1][3])) + Q[0][1];
	Pcap[0][2] = ((FP[0][0]*F[2][0])+(FP[0][1]*F[2][1])+(FP[0][2]*F[2][2])+(FP[0][3]*F[2][3])) + Q[0][2];
	Pcap[0][3] = ((FP[0][0]*F[3][0])+(FP[0][1]*F[3][1])+(FP[0][2]*F[3][2])+(FP[0][3]*F[3][3])) + Q[0][3];

	Pcap[1][0] = ((FP[1][0]*F[0][0])+(FP[1][1]*F[0][1])+(FP[1][2]*F[0][2])+(FP[1][3]*F[0][3])) + Q[1][0];
	Pcap[1][1] = ((FP[1][0]*F[1][0])+(FP[1][1]*F[1][1])+(FP[1][2]*F[1][2])+(FP[1][3]*F[1][3])) + Q[1][1];
	Pcap[1][2] = ((FP[1][0]*F[2][0])+(FP[1][1]*F[2][1])+(FP[1][2]*F[2][2])+(FP[1][3]*F[2][3])) + Q[1][2];
	Pcap[1][3] = ((FP[1][0]*F[3][0])+(FP[1][1]*F[3][1])+(FP[1][2]*F[3][2])+(FP[1][3]*F[3][3])) + Q[1][3];

	Pcap[2][0] = ((FP[2][0]*F[0][0])+(FP[2][1]*F[0][1])+(FP[2][2]*F[0][2])+(FP[2][3]*F[0][3])) + Q[2][0];
	Pcap[2][1] = ((FP[2][0]*F[1][0])+(FP[2][1]*F[1][1])+(FP[2][2]*F[1][2])+(FP[2][3]*F[1][3])) + Q[2][1];
	Pcap[2][2] = ((FP[2][0]*F[2][0])+(FP[2][1]*F[2][1])+(FP[2][2]*F[2][2])+(FP[2][3]*F[2][3])) + Q[2][2];
	Pcap[2][3] = ((FP[2][0]*F[3][0])+(FP[2][1]*F[3][1])+(FP[2][2]*F[3][2])+(FP[2][3]*F[3][3])) + Q[2][3];

	Pcap[3][0] = ((FP[3][0]*F[0][0])+(FP[3][1]*F[0][1])+(FP[3][2]*F[0][2])+(FP[3][3]*F[0][3])) + Q[3][0];
	Pcap[3][1] = ((FP[3][0]*F[1][0])+(FP[3][1]*F[1][1])+(FP[3][2]*F[1][2])+(FP[3][3]*F[1][3])) + Q[3][1];
	Pcap[3][2] = ((FP[3][0]*F[2][0])+(FP[3][1]*F[2][1])+(FP[3][2]*F[2][2])+(FP[3][3]*F[2][3])) + Q[3][2];
	Pcap[3][3] = ((FP[3][0]*F[3][0])+(FP[3][1]*F[3][1])+(FP[3][2]*F[3][2])+(FP[3][3]*F[3][3])) + Q[3][3];

	// Correction: q = qcap + Kg(z - h(q))

	// acap = Cqt * g (NED frame)
	acap[0] = 2*(qcap.x*qcap.z-qcap.s*qcap.y);
	acap[1] = 2*(qcap.s*qcap.x+qcap.y*qcap.z);
	acap[2] = qcapS_Squared-qcapX_Squared-qcapY_Squared+qcapZ_Squared;


	// mcap = Cqt * r (NED frame)
	mcap[0] = qcapS_Squared+qcapX_Squared-qcapY_Squared-qcapZ_Squared*rx + 2*(qcap.x*qcap.z-qcap.s*qcap.y)*rz;
	mcap[1] = 2*(qcap.x*qcap.y-qcap.s*qcap.z)*rx + 2*(qcap.s*qcap.x+qcap.y*qcap.z)*rz;
	mcap[2] = 2*(qcap.x*qcap.z+qcap.s*qcap.y)*rx + qcapS_Squared-qcapX_Squared-qcapY_Squared+qcapZ_Squared*rz;

	// v = z - h(qcap)
	// z = [a, m]^T and h(qcap) = [acap, mcap]^T
	v[0] = ax - acap[0];
	v[1] = ay - acap[1];
	v[2] = az - acap[2];
	v[3] = mx - mcap[0];
	v[4] = my - mcap[1];
	v[5] = mz - mcap[2];

	// S = H(qcap)*Pcap*H(qcap)^T + R
	// H(qcap) = Jacobian of h(qcap)
	Hqcap[0][0] = 2 * -qcap.y;	Hqcap[0][1] = 2 * qcap.z;	Hqcap[0][2] = 2 * -qcap.s;	Hqcap[0][2] = 2 * -qcap.x;
	Hqcap[1][0] = 2 * qcap.x;	Hqcap[1][1] = 2 * qcap.s;	Hqcap[1][2] = 2 * qcap.z;	Hqcap[1][2] = 2 * qcap.y;
	Hqcap[2][0] = 2 * qcap.s;	Hqcap[2][1] = 2 * -qcap.x;	Hqcap[2][2] = 2 * -qcap.y;	Hqcap[2][2] = 2 * qcap.z;
	Hqcap[3][0] = 2 * (rx*q.s - rz*q.y);  Hqcap[3][1] = 2 * (rx*q.x + rz*q.z);  Hqcap[3][2] = 2 * (-rx*q.y - rz*q.s); Hqcap[3][3] = 2 * (-rx*q.z + rz*q.x);
	Hqcap[4][0] = 2 * (-rx*q.z + rz*q.x); Hqcap[4][1] = 2 * (rx*q.y + rz*q.s);  Hqcap[4][2] = 2 * (rx*q.s - rz*q.y);  Hqcap[4][3] = 2 * (-rx*q.x - rz*q.z);
	Hqcap[5][0] = 2 * (rx*q.y - rz*q.s);  Hqcap[5][1] = 2 * (rx*q.z - rz*q.x);  Hqcap[5][2] = 2 * (rx*q.x + rz*q.z);  Hqcap[5][3] = 2 * (rx*q.s - rz*q.y);

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

//	HqcapPcap[0][0] = Hqcap[0][0]*Pcap[0][0]+Hqcap[0][1]*Pcap[1][0]+Hqcap[0][2]*Pcap[2][0]+Hqcap[0][3]*Pcap[3][0];
//	HqcapPcap[0][1] = Hqcap[0][0]*Pcap[0][1]+Hqcap[0][1]*Pcap[1][1]+Hqcap[0][2]*Pcap[2][1]+Hqcap[0][3]*Pcap[3][1];
//	HqcapPcap[0][2] = Hqcap[0][0]*Pcap[0][2]+Hqcap[0][1]*Pcap[1][2]+Hqcap[0][2]*Pcap[2][2]+Hqcap[0][3]*Pcap[3][2];
//	HqcapPcap[0][3] = Hqcap[0][0]*Pcap[0][3]+Hqcap[0][1]*Pcap[1][3]+Hqcap[0][2]*Pcap[2][3]+Hqcap[0][3]*Pcap[3][3];
//	HqcapPcap[1][0] = Hqcap[1][0]*Pcap[0][0]+Hqcap[1][1]*Pcap[1][0]+Hqcap[1][2]*Pcap[2][0]+Hqcap[1][3]*Pcap[3][0];
//	HqcapPcap[1][1] = Hqcap[1][0]*Pcap[0][1]+Hqcap[1][1]*Pcap[1][1]+Hqcap[1][2]*Pcap[2][1]+Hqcap[1][3]*Pcap[3][1];
//	HqcapPcap[1][2] = Hqcap[1][0]*Pcap[0][2]+Hqcap[1][1]*Pcap[1][2]+Hqcap[1][2]*Pcap[2][2]+Hqcap[1][3]*Pcap[3][2];
//	HqcapPcap[1][3] = Hqcap[1][0]*Pcap[0][3]+Hqcap[1][1]*Pcap[1][3]+Hqcap[1][2]*Pcap[2][3]+Hqcap[1][3]*Pcap[3][3];
//	HqcapPcap[2][0] = Hqcap[2][0]*Pcap[0][0]+Hqcap[2][1]*Pcap[1][0]+Hqcap[2][2]*Pcap[2][0]+Hqcap[2][3]*Pcap[3][0];
//	HqcapPcap[2][1] = Hqcap[2][0]*Pcap[0][1]+Hqcap[2][1]*Pcap[1][1]+Hqcap[2][2]*Pcap[2][1]+Hqcap[2][3]*Pcap[3][1];
//	HqcapPcap[2][2] = Hqcap[2][0]*Pcap[0][2]+Hqcap[2][1]*Pcap[1][2]+Hqcap[2][2]*Pcap[2][2]+Hqcap[2][3]*Pcap[3][2];
//	HqcapPcap[2][3] = Hqcap[2][0]*Pcap[0][3]+Hqcap[2][1]*Pcap[1][3]+Hqcap[2][2]*Pcap[2][3]+Hqcap[2][3]*Pcap[3][3];
//	HqcapPcap[3][0] = Hqcap[3][0]*Pcap[0][0]+Hqcap[3][1]*Pcap[1][0]+Hqcap[3][2]*Pcap[2][0]+Hqcap[3][3]*Pcap[3][0];
//	HqcapPcap[3][1] = Hqcap[3][0]*Pcap[0][1]+Hqcap[3][1]*Pcap[1][1]+Hqcap[3][2]*Pcap[2][1]+Hqcap[3][3]*Pcap[3][1];
//	HqcapPcap[3][2] = Hqcap[3][0]*Pcap[0][2]+Hqcap[3][1]*Pcap[1][2]+Hqcap[3][2]*Pcap[2][2]+Hqcap[3][3]*Pcap[3][2];
//	HqcapPcap[3][3] = Hqcap[3][0]*Pcap[0][3]+Hqcap[3][1]*Pcap[1][3]+Hqcap[3][2]*Pcap[2][3]+Hqcap[3][3]*Pcap[3][3];
//	HqcapPcap[4][0] = Hqcap[4][0]*Pcap[0][0]+Hqcap[4][1]*Pcap[1][0]+Hqcap[4][2]*Pcap[2][0]+Hqcap[4][3]*Pcap[3][0];
//	HqcapPcap[4][1] = Hqcap[4][0]*Pcap[0][1]+Hqcap[4][1]*Pcap[1][1]+Hqcap[4][2]*Pcap[2][1]+Hqcap[4][3]*Pcap[3][1];
//	HqcapPcap[4][2] = Hqcap[4][0]*Pcap[0][2]+Hqcap[4][1]*Pcap[1][2]+Hqcap[4][2]*Pcap[2][2]+Hqcap[4][3]*Pcap[3][2];
//	HqcapPcap[4][3] = Hqcap[4][0]*Pcap[0][3]+Hqcap[4][1]*Pcap[1][3]+Hqcap[4][2]*Pcap[2][3]+Hqcap[4][3]*Pcap[3][3];
//	HqcapPcap[5][0] = Hqcap[5][0]*Pcap[0][0]+Hqcap[5][1]*Pcap[1][0]+Hqcap[5][2]*Pcap[2][0]+Hqcap[5][3]*Pcap[3][0];
//	HqcapPcap[5][1] = Hqcap[5][0]*Pcap[0][1]+Hqcap[5][1]*Pcap[1][1]+Hqcap[5][2]*Pcap[2][1]+Hqcap[5][3]*Pcap[3][1];
//	HqcapPcap[5][2] = Hqcap[5][0]*Pcap[0][2]+Hqcap[5][1]*Pcap[1][2]+Hqcap[5][2]*Pcap[2][2]+Hqcap[5][3]*Pcap[3][2];
//	HqcapPcap[5][3] = Hqcap[5][0]*Pcap[0][3]+Hqcap[5][1]*Pcap[1][3]+Hqcap[5][2]*Pcap[2][3]+Hqcap[5][3]*Pcap[3][3];

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

//	S[0][0] = HqcapPcap[0][0]*Hqcap[0][0]+HqcapPcap[0][1]*Hqcap[0][1]+HqcapPcap[0][2]*Hqcap[0][2]+HqcapPcap[0][3]*Hqcap[0][3];
//	S[0][1] = HqcapPcap[0][0]*Hqcap[1][0]+HqcapPcap[0][1]*Hqcap[1][1]+HqcapPcap[0][2]*Hqcap[1][2]+HqcapPcap[0][3]*Hqcap[1][3];
//	S[0][2] = HqcapPcap[0][0]*Hqcap[2][0]+HqcapPcap[0][1]*Hqcap[2][1]+HqcapPcap[0][2]*Hqcap[2][2]+HqcapPcap[0][3]*Hqcap[2][3];
//	S[0][3] = HqcapPcap[0][0]*Hqcap[3][0]+HqcapPcap[0][1]*Hqcap[3][1]+HqcapPcap[0][2]*Hqcap[3][2]+HqcapPcap[0][3]*Hqcap[3][3];
//	S[0][4] = HqcapPcap[0][0]*Hqcap[4][0]+HqcapPcap[0][1]*Hqcap[4][1]+HqcapPcap[0][2]*Hqcap[4][2]+HqcapPcap[0][3]*Hqcap[4][3];
//	S[0][5] = HqcapPcap[0][0]*Hqcap[5][0]+HqcapPcap[0][1]*Hqcap[5][1]+HqcapPcap[0][2]*Hqcap[5][2]+HqcapPcap[0][3]*Hqcap[5][3];
//	S[1][0] = HqcapPcap[1][0]*Hqcap[0][0]+HqcapPcap[1][1]*Hqcap[0][1]+HqcapPcap[1][2]*Hqcap[0][2]+HqcapPcap[1][3]*Hqcap[0][3];
//	S[1][1] = HqcapPcap[1][0]*Hqcap[1][0]+HqcapPcap[1][1]*Hqcap[1][1]+HqcapPcap[1][2]*Hqcap[1][2]+HqcapPcap[1][3]*Hqcap[1][3];
//	S[1][2] = HqcapPcap[1][0]*Hqcap[2][0]+HqcapPcap[1][1]*Hqcap[2][1]+HqcapPcap[1][2]*Hqcap[2][2]+HqcapPcap[1][3]*Hqcap[2][3];
//	S[1][3] = HqcapPcap[1][0]*Hqcap[3][0]+HqcapPcap[1][1]*Hqcap[3][1]+HqcapPcap[1][2]*Hqcap[3][2]+HqcapPcap[1][3]*Hqcap[3][3];
//	S[1][4] = HqcapPcap[1][0]*Hqcap[4][0]+HqcapPcap[1][1]*Hqcap[4][1]+HqcapPcap[1][2]*Hqcap[4][2]+HqcapPcap[1][3]*Hqcap[4][3];
//	S[1][5] = HqcapPcap[1][0]*Hqcap[5][0]+HqcapPcap[1][1]*Hqcap[5][1]+HqcapPcap[1][2]*Hqcap[5][2]+HqcapPcap[1][3]*Hqcap[5][3];
//	S[2][0] = HqcapPcap[2][0]*Hqcap[0][0]+HqcapPcap[2][1]*Hqcap[0][1]+HqcapPcap[2][2]*Hqcap[0][2]+HqcapPcap[2][3]*Hqcap[0][3];
//	S[2][1] = HqcapPcap[2][0]*Hqcap[1][0]+HqcapPcap[2][1]*Hqcap[1][1]+HqcapPcap[2][2]*Hqcap[1][2]+HqcapPcap[2][3]*Hqcap[1][3];
//	S[2][2] = HqcapPcap[2][0]*Hqcap[2][0]+HqcapPcap[2][1]*Hqcap[2][1]+HqcapPcap[2][2]*Hqcap[2][2]+HqcapPcap[2][3]*Hqcap[2][3];
//	S[2][3] = HqcapPcap[2][0]*Hqcap[3][0]+HqcapPcap[2][1]*Hqcap[3][1]+HqcapPcap[2][2]*Hqcap[3][2]+HqcapPcap[2][3]*Hqcap[3][3];
//	S[2][4] = HqcapPcap[2][0]*Hqcap[4][0]+HqcapPcap[2][1]*Hqcap[4][1]+HqcapPcap[2][2]*Hqcap[4][2]+HqcapPcap[2][3]*Hqcap[4][3];
//	S[2][5] = HqcapPcap[2][0]*Hqcap[5][0]+HqcapPcap[2][1]*Hqcap[5][1]+HqcapPcap[2][2]*Hqcap[5][2]+HqcapPcap[2][3]*Hqcap[5][3];
//	S[3][0] = HqcapPcap[3][0]*Hqcap[0][0]+HqcapPcap[3][1]*Hqcap[0][1]+HqcapPcap[3][2]*Hqcap[0][2]+HqcapPcap[3][3]*Hqcap[0][3];
//	S[3][1] = HqcapPcap[3][0]*Hqcap[1][0]+HqcapPcap[3][1]*Hqcap[1][1]+HqcapPcap[3][2]*Hqcap[1][2]+HqcapPcap[3][3]*Hqcap[1][3];
//	S[3][2] = HqcapPcap[3][0]*Hqcap[2][0]+HqcapPcap[3][1]*Hqcap[2][1]+HqcapPcap[3][2]*Hqcap[2][2]+HqcapPcap[3][3]*Hqcap[2][3];
//	S[3][3] = HqcapPcap[3][0]*Hqcap[3][0]+HqcapPcap[3][1]*Hqcap[3][1]+HqcapPcap[3][2]*Hqcap[3][2]+HqcapPcap[3][3]*Hqcap[3][3];
//	S[3][4] = HqcapPcap[3][0]*Hqcap[4][0]+HqcapPcap[3][1]*Hqcap[4][1]+HqcapPcap[3][2]*Hqcap[4][2]+HqcapPcap[3][3]*Hqcap[4][3];
//	S[3][5] = HqcapPcap[3][0]*Hqcap[5][0]+HqcapPcap[3][1]*Hqcap[5][1]+HqcapPcap[3][2]*Hqcap[5][2]+HqcapPcap[3][3]*Hqcap[5][3];
//	S[4][0] = HqcapPcap[4][0]*Hqcap[0][0]+HqcapPcap[4][1]*Hqcap[0][1]+HqcapPcap[4][2]*Hqcap[0][2]+HqcapPcap[4][3]*Hqcap[0][3];
//	S[4][1] = HqcapPcap[4][0]*Hqcap[1][0]+HqcapPcap[4][1]*Hqcap[1][1]+HqcapPcap[4][2]*Hqcap[1][2]+HqcapPcap[4][3]*Hqcap[1][3];
//	S[4][2] = HqcapPcap[4][0]*Hqcap[2][0]+HqcapPcap[4][1]*Hqcap[2][1]+HqcapPcap[4][2]*Hqcap[2][2]+HqcapPcap[4][3]*Hqcap[2][3];
//	S[4][3] = HqcapPcap[4][0]*Hqcap[3][0]+HqcapPcap[4][1]*Hqcap[3][1]+HqcapPcap[4][2]*Hqcap[3][2]+HqcapPcap[4][3]*Hqcap[3][3];
//	S[4][4] = HqcapPcap[4][0]*Hqcap[4][0]+HqcapPcap[4][1]*Hqcap[4][1]+HqcapPcap[4][2]*Hqcap[4][2]+HqcapPcap[4][3]*Hqcap[4][3];
//	S[4][5] = HqcapPcap[4][0]*Hqcap[5][0]+HqcapPcap[4][1]*Hqcap[5][1]+HqcapPcap[4][2]*Hqcap[5][2]+HqcapPcap[4][3]*Hqcap[5][3];
//	S[5][0] = HqcapPcap[5][0]*Hqcap[0][0]+HqcapPcap[5][1]*Hqcap[0][1]+HqcapPcap[5][2]*Hqcap[0][2]+HqcapPcap[5][3]*Hqcap[0][3];
//	S[5][1] = HqcapPcap[5][0]*Hqcap[1][0]+HqcapPcap[5][1]*Hqcap[1][1]+HqcapPcap[5][2]*Hqcap[1][2]+HqcapPcap[5][3]*Hqcap[1][3];
//	S[5][2] = HqcapPcap[5][0]*Hqcap[2][0]+HqcapPcap[5][1]*Hqcap[2][1]+HqcapPcap[5][2]*Hqcap[2][2]+HqcapPcap[5][3]*Hqcap[2][3];
//	S[5][3] = HqcapPcap[5][0]*Hqcap[3][0]+HqcapPcap[5][1]*Hqcap[3][1]+HqcapPcap[5][2]*Hqcap[3][2]+HqcapPcap[5][3]*Hqcap[3][3];
//	S[5][4] = HqcapPcap[5][0]*Hqcap[4][0]+HqcapPcap[5][1]*Hqcap[4][1]+HqcapPcap[5][2]*Hqcap[4][2]+HqcapPcap[5][3]*Hqcap[4][3];
//	S[5][5] = HqcapPcap[5][0]*Hqcap[5][0]+HqcapPcap[5][1]*Hqcap[5][1]+HqcapPcap[5][2]*Hqcap[5][2]+HqcapPcap[5][3]*Hqcap[5][3];


	// H(qcap)*Pcap*H(qcap)^T + R = S + R
	// R is measurement noise of the accelerometer and magnetometer, it is a 6x6 diagonal matrix
	S[0][0] += R[0];
	S[1][1] += R[1];
	S[2][2] += R[2];
	S[3][3] += R[3];
	S[4][4] += R[4];
	S[5][5] += R[5];

	// calculating inverse of S matrix, using LU decomposition
	// S * S^-1 = I, using LU decomposition S = LU
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
 			for (uint8_t c = 0; c < 6; c++)
 				sum += s[i][c] * Y[c];
 		}

 		// Since, U * X = Y, here X = one column of S^-1
        // Now, solve for S^-1 (sinv)
 		for (uint8_t i = 5; i > -1; i--)
 		{
 			sinv[i][col] = Y[i];
 			for (uint8_t j = i + 1; j < 6; j++)
 				sinv[i][col] -= s[i][j] * sinv[j][col];

 			// singular pivot
 			if (fabs(s[i][i]) <= 0.00001f)
 				return false;

 			sinv[i][col] /= s[i][i];
 		}
	}

 	// Kg = Pcap * HT * S^-1(sinv)
 	// here, Kg = PcapHT * sinv
 	Kg[0][0] = PcapHT[0][0]*sinv[0][0]+PcapHT[0][1]*sinv[1][0]+PcapHT[0][2]*sinv[2][0]+PcapHT[0][3]*sinv[3][0];
 	Kg[0][1] = PcapHT[0][0]*sinv[0][1]+PcapHT[0][1]*sinv[1][1]+PcapHT[0][2]*sinv[2][1]+PcapHT[0][3]*sinv[3][1];
 	Kg[0][2] = PcapHT[0][0]*sinv[0][2]+PcapHT[0][1]*sinv[1][2]+PcapHT[0][2]*sinv[2][2]+PcapHT[0][3]*sinv[3][2];
 	Kg[0][3] = PcapHT[0][0]*sinv[0][3]+PcapHT[0][1]*sinv[1][3]+PcapHT[0][2]*sinv[2][3]+PcapHT[0][3]*sinv[3][3];
 	Kg[0][4] = PcapHT[0][0]*sinv[0][4]+PcapHT[0][1]*sinv[1][4]+PcapHT[0][2]*sinv[2][4]+PcapHT[0][3]*sinv[3][4];
 	Kg[0][5] = PcapHT[0][0]*sinv[0][5]+PcapHT[0][1]*sinv[1][5]+PcapHT[0][2]*sinv[2][5]+PcapHT[0][3]*sinv[3][5];
 	Kg[1][0] = PcapHT[1][0]*sinv[0][0]+PcapHT[1][1]*sinv[1][0]+PcapHT[1][2]*sinv[2][0]+PcapHT[1][3]*sinv[3][0];
 	Kg[1][1] = PcapHT[1][0]*sinv[0][1]+PcapHT[1][1]*sinv[1][1]+PcapHT[1][2]*sinv[2][1]+PcapHT[1][3]*sinv[3][1];
 	Kg[1][2] = PcapHT[1][0]*sinv[0][2]+PcapHT[1][1]*sinv[1][2]+PcapHT[1][2]*sinv[2][2]+PcapHT[1][3]*sinv[3][2];
 	Kg[1][3] = PcapHT[1][0]*sinv[0][3]+PcapHT[1][1]*sinv[1][3]+PcapHT[1][2]*sinv[2][3]+PcapHT[1][3]*sinv[3][3];
 	Kg[1][4] = PcapHT[1][0]*sinv[0][4]+PcapHT[1][1]*sinv[1][4]+PcapHT[1][2]*sinv[2][4]+PcapHT[1][3]*sinv[3][4];
 	Kg[1][5] = PcapHT[1][0]*sinv[0][5]+PcapHT[1][1]*sinv[1][5]+PcapHT[1][2]*sinv[2][5]+PcapHT[1][3]*sinv[3][5];
 	Kg[2][0] = PcapHT[2][0]*sinv[0][0]+PcapHT[2][1]*sinv[1][0]+PcapHT[2][2]*sinv[2][0]+PcapHT[2][3]*sinv[3][0];
 	Kg[2][1] = PcapHT[2][0]*sinv[0][1]+PcapHT[2][1]*sinv[1][1]+PcapHT[2][2]*sinv[2][1]+PcapHT[2][3]*sinv[3][1];
 	Kg[2][2] = PcapHT[2][0]*sinv[0][2]+PcapHT[2][1]*sinv[1][2]+PcapHT[2][2]*sinv[2][2]+PcapHT[2][3]*sinv[3][2];
 	Kg[2][3] = PcapHT[2][0]*sinv[0][3]+PcapHT[2][1]*sinv[1][3]+PcapHT[2][2]*sinv[2][3]+PcapHT[2][3]*sinv[3][3];
 	Kg[2][4] = PcapHT[2][0]*sinv[0][4]+PcapHT[2][1]*sinv[1][4]+PcapHT[2][2]*sinv[2][4]+PcapHT[2][3]*sinv[3][4];
 	Kg[2][5] = PcapHT[2][0]*sinv[0][5]+PcapHT[2][1]*sinv[1][5]+PcapHT[2][2]*sinv[2][5]+PcapHT[2][3]*sinv[3][5];
 	Kg[3][0] = PcapHT[3][0]*sinv[0][0]+PcapHT[3][1]*sinv[1][0]+PcapHT[3][2]*sinv[2][0]+PcapHT[3][3]*sinv[3][0];
 	Kg[3][1] = PcapHT[3][0]*sinv[0][1]+PcapHT[3][1]*sinv[1][1]+PcapHT[3][2]*sinv[2][1]+PcapHT[3][3]*sinv[3][1];
 	Kg[3][2] = PcapHT[3][0]*sinv[0][2]+PcapHT[3][1]*sinv[1][2]+PcapHT[3][2]*sinv[2][2]+PcapHT[3][3]*sinv[3][2];
 	Kg[3][3] = PcapHT[3][0]*sinv[0][3]+PcapHT[3][1]*sinv[1][3]+PcapHT[3][2]*sinv[2][3]+PcapHT[3][3]*sinv[3][3];
 	Kg[3][4] = PcapHT[3][0]*sinv[0][4]+PcapHT[3][1]*sinv[1][4]+PcapHT[3][2]*sinv[2][4]+PcapHT[3][3]*sinv[3][4];
 	Kg[3][5] = PcapHT[3][0]*sinv[0][5]+PcapHT[3][1]*sinv[1][5]+PcapHT[3][2]*sinv[2][5]+PcapHT[3][3]*sinv[3][5];

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

 	// P = (I4 - Kg * H(qcap)) * Pcap
 	// Solving P = (I4 - Kg * Hqcap)
 	P[0][0] = 1.0 - (Kg[0][0]*Hqcap[0][0]+Kg[0][1]*Hqcap[1][0]+Kg[0][2]*Hqcap[2][0]+Kg[0][3]*Hqcap[3][0]+Kg[0][4]*Hqcap[4][0]+Kg[0][5]*Hqcap[5][0]);
 	P[0][1] = -(Kg[0][0]*Hqcap[0][1]+Kg[0][1]*Hqcap[1][1]+Kg[0][2]*Hqcap[2][1]+Kg[0][3]*Hqcap[3][1]+Kg[0][4]*Hqcap[4][1]+Kg[0][5]*Hqcap[5][1]);
 	P[0][2] = -(Kg[0][0]*Hqcap[0][2]+Kg[0][1]*Hqcap[1][2]+Kg[0][2]*Hqcap[2][2]+Kg[0][3]*Hqcap[3][2]+Kg[0][4]*Hqcap[4][2]+Kg[0][5]*Hqcap[5][2]);
 	P[0][3] = -(Kg[0][0]*Hqcap[0][3]+Kg[0][1]*Hqcap[1][3]+Kg[0][2]*Hqcap[2][3]+Kg[0][3]*Hqcap[3][3]+Kg[0][4]*Hqcap[4][3]+Kg[0][5]*Hqcap[5][3]);
 	P[1][0] = -(Kg[1][0]*Hqcap[0][0]+Kg[1][1]*Hqcap[1][0]+Kg[1][2]*Hqcap[2][0]+Kg[1][3]*Hqcap[3][0]+Kg[1][4]*Hqcap[4][0]+Kg[1][5]*Hqcap[5][0]);
 	P[1][1] = 1.0 - (Kg[1][0]*Hqcap[0][1]+Kg[1][1]*Hqcap[1][1]+Kg[1][2]*Hqcap[2][1]+Kg[1][3]*Hqcap[3][1]+Kg[1][4]*Hqcap[4][1]+Kg[1][5]*Hqcap[5][1]);
 	P[1][2] = -(Kg[1][0]*Hqcap[0][2]+Kg[1][1]*Hqcap[1][2]+Kg[1][2]*Hqcap[2][2]+Kg[1][3]*Hqcap[3][2]+Kg[1][4]*Hqcap[4][2]+Kg[1][5]*Hqcap[5][2]);
 	P[1][3] = -(Kg[1][0]*Hqcap[0][3]+Kg[1][1]*Hqcap[1][3]+Kg[1][2]*Hqcap[2][3]+Kg[1][3]*Hqcap[3][3]+Kg[1][4]*Hqcap[4][3]+Kg[1][5]*Hqcap[5][3]);
 	P[2][0] = -(Kg[2][0]*Hqcap[0][0]+Kg[2][1]*Hqcap[1][0]+Kg[2][2]*Hqcap[2][0]+Kg[2][3]*Hqcap[3][0]+Kg[2][4]*Hqcap[4][0]+Kg[2][5]*Hqcap[5][0]);
 	P[2][1] = -(Kg[2][0]*Hqcap[0][1]+Kg[2][1]*Hqcap[1][1]+Kg[2][2]*Hqcap[2][1]+Kg[2][3]*Hqcap[3][1]+Kg[2][4]*Hqcap[4][1]+Kg[2][5]*Hqcap[5][1]);
 	P[2][2] = 1.0 - (Kg[2][0]*Hqcap[0][2]+Kg[2][1]*Hqcap[1][2]+Kg[2][2]*Hqcap[2][2]+Kg[2][3]*Hqcap[3][2]+Kg[2][4]*Hqcap[4][2]+Kg[2][5]*Hqcap[5][2]);
 	P[2][3] = -(Kg[2][0]*Hqcap[0][3]+Kg[2][1]*Hqcap[1][3]+Kg[2][2]*Hqcap[2][3]+Kg[2][3]*Hqcap[3][3]+Kg[2][4]*Hqcap[4][3]+Kg[2][5]*Hqcap[5][3]);
 	P[3][0] = -(Kg[3][0]*Hqcap[0][0]+Kg[3][1]*Hqcap[1][0]+Kg[3][2]*Hqcap[2][0]+Kg[3][3]*Hqcap[3][0]+Kg[3][4]*Hqcap[4][0]+Kg[3][5]*Hqcap[5][0]);
 	P[3][1] = -(Kg[3][0]*Hqcap[0][1]+Kg[3][1]*Hqcap[1][1]+Kg[3][2]*Hqcap[2][1]+Kg[3][3]*Hqcap[3][1]+Kg[3][4]*Hqcap[4][1]+Kg[3][5]*Hqcap[5][1]);
 	P[3][2] = -(Kg[3][0]*Hqcap[0][2]+Kg[3][1]*Hqcap[1][2]+Kg[3][2]*Hqcap[2][2]+Kg[3][3]*Hqcap[3][2]+Kg[3][4]*Hqcap[4][2]+Kg[3][5]*Hqcap[5][2]);
 	P[3][3] = 1.0 - (Kg[3][0]*Hqcap[0][3]+Kg[3][1]*Hqcap[1][3]+Kg[3][2]*Hqcap[2][3]+Kg[3][3]*Hqcap[3][3]+Kg[3][4]*Hqcap[4][3]+Kg[3][5]*Hqcap[5][3]);

 	// Now, P = P * Pcap
 	P[0][0] = P[0][0]*Pcap[0][0]+P[0][1]*Pcap[1][0]+P[0][2]*Pcap[2][0]+P[0][3]*Pcap[3][0];
 	P[0][1] = P[0][0]*Pcap[0][1]+P[0][1]*Pcap[1][1]+P[0][2]*Pcap[2][1]+P[0][3]*Pcap[3][1];
 	P[0][2] = P[0][0]*Pcap[0][2]+P[0][1]*Pcap[1][2]+P[0][2]*Pcap[2][2]+P[0][3]*Pcap[3][2];
 	P[0][3] = P[0][0]*Pcap[0][3]+P[0][1]*Pcap[1][3]+P[0][2]*Pcap[2][3]+P[0][3]*Pcap[3][3];
 	P[1][0] = P[1][0]*Pcap[0][0]+P[1][1]*Pcap[1][0]+P[1][2]*Pcap[2][0]+P[1][3]*Pcap[3][0];
 	P[1][1] = P[1][0]*Pcap[0][1]+P[1][1]*Pcap[1][1]+P[1][2]*Pcap[2][1]+P[1][3]*Pcap[3][1];
 	P[1][2] = P[1][0]*Pcap[0][2]+P[1][1]*Pcap[1][2]+P[1][2]*Pcap[2][2]+P[1][3]*Pcap[3][2];
 	P[1][3] = P[1][0]*Pcap[0][3]+P[1][1]*Pcap[1][3]+P[1][2]*Pcap[2][3]+P[1][3]*Pcap[3][3];
 	P[2][0] = P[2][0]*Pcap[0][0]+P[2][1]*Pcap[1][0]+P[2][2]*Pcap[2][0]+P[2][3]*Pcap[3][0];
 	P[2][1] = P[2][0]*Pcap[0][1]+P[2][1]*Pcap[1][1]+P[2][2]*Pcap[2][1]+P[2][3]*Pcap[3][1];
 	P[2][2] = P[2][0]*Pcap[0][2]+P[2][1]*Pcap[1][2]+P[2][2]*Pcap[2][2]+P[2][3]*Pcap[3][2];
 	P[2][3] = P[2][0]*Pcap[0][3]+P[2][1]*Pcap[1][3]+P[2][2]*Pcap[2][3]+P[2][3]*Pcap[3][3];
 	P[3][0] = P[3][0]*Pcap[0][0]+P[3][1]*Pcap[1][0]+P[3][2]*Pcap[2][0]+P[3][3]*Pcap[3][0];
 	P[3][1] = P[3][0]*Pcap[0][1]+P[3][1]*Pcap[1][1]+P[3][2]*Pcap[2][1]+P[3][3]*Pcap[3][1];
 	P[3][2] = P[3][0]*Pcap[0][2]+P[3][1]*Pcap[1][2]+P[3][2]*Pcap[2][2]+P[3][3]*Pcap[3][2];
 	P[3][3] = P[3][0]*Pcap[0][3]+P[3][1]*Pcap[1][3]+P[3][2]*Pcap[2][3]+P[3][3]*Pcap[3][3];

 	q.Normalise();

 	return true;
}
