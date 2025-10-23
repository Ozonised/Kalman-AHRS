#ifndef EXTENDEDKALMANFILTER_HPP_
#define _EXTENDEDKALMANFILTER_HPP_

#include "quaternion.hpp"
#include "math.h"
#include <stdbool.h>
#include <string.h>
#include <stdint.h>

class ExtendedKalmanFilter
{
	private:
		Quaternion qcap; 		// next state estimate
		Quaternion q;			// state / output
		float v[6];				// Measurement Residual/ Error / Difference between measurement(z) and estimate(Xcap)
		float acap[3];
		float mcap[3];
		float Kg[4][6];			// Kalman Gain
		float Pcap[4][4];		// Predicted Error State Covariance (before measurement z)
		float P[4][4];			// Estimated Covariance of the state (after measurement z)
		float F[4][4];			// Fundamental Matrix/ State Transition Matrix
		float Hqcap[6][4];
		float PcapHT[4][6];
		float S[6][6];			// Measurement Prediction Covariance (also stores S^-1)
		float s[6][6];			// copy of S
		float SigmaOmega[3]; 	// Gyro spectral noise covariance
		float R[6];				// Measurement noise covariance matrix
		float Q[4][4];
		float dt;				// sampling time
		float dt2;				// half of sampling time
		float rx, rz;			// magnetic inclination
		float magDeclination;	// magnetic declination

	public:
		ExtendedKalmanFilter();
		void SetSampleTime(float freq);
		void SetInitialState(float ax, float ay, float az, float mx, float my, float mz);
		void SetMagneticDip(float degrees);
		void SetGyroNoise(float NoiseX, float NoiseY, float NoiseZ);
		void SetR(float NoiseAx, float NoiseAy, float NoiseAz, float NoiseMx, float NoiseMy, float NoiseMz);
		bool Run(float ax, float ay, float az, float gx, float gy, float gz, float mx, float my, float mz);
		void GetOrientation(Quaternion& qState);
};

#endif /*EXTENDEDKALMANFILTER_HPP_ */
