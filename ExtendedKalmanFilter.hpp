#ifndef EXTENDEDKALMANFILTER_HPP_
#define _EXTENDEDKALMANFILTER_HPP_

#include "quanternion.hpp"
#include "math.h"
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
// default magnetometer noise
class ExtendedKalmanFilter
{
	private:
		Quanternion qcap; 		// next state estimate
		Quanternion q;			// state / output
		float v[6];				// Measurement Residual/ Error / Difference between measurement(z) and estimate(Xcap)
		float U[3];				// gyroscope readings
		float a[3];
		float m[3];
		float acap[3];
		float mcap[3];
		float Kg[4][6];			// Kalman Gain
		float Pcap[4][4];		// Predicted Error State Covariance (before measurement z)
		float P[4][4];			// Estimated Covariance of the state (after measurement z)
		float F[4][4];			// Fundamental Matrix/ State Transition Matrix
		float FP[4][4];			// Product of F and P
		float Q[4][4];			// Process Noise Covariance Matrix
		float Hqcap[6][4];
		float PcapHT[4][6];
		float S[6][6];			// Measurement Prediction Covariance
		float s[6][6];			// copy of S
		float sinv[6][6];		// S^-1
		float Icol[6], Y[6];
		float SigmaOmega[3]; 	// Gyro spectral noise covariance
		float R[6];				// Measurement noise covariance matrix
		float dt2;				// half of sampling time
		float rx, rz;			// magnetic inclination

	public:
		ExtendedKalmanFilter();
		void UpdateU(float gx, float gy, float gz);
		void SetGyroNoise(float x, float y, float z);
		void SetR(float NoiseAx, float NoiseAy, float NoiseAz, float NoiseMx, float NoiseMy, float NoiseMz);
		void SetMagneticInclination(float degrees);
		bool Run(float ax, float ay, float az, float mz, float gx, float gy, float gz, float mx, float my);
		void SetSampleTime(float dt2);

};

#endif /*EXTENDEDKALMANFILTER_HPP_ */
