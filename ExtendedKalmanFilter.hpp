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
		Quaternion qam;
		float Kg[4][4];			// Kalman Gain
		float Pcap[4][4];		// Predicted Error State Covariance (before measurement z)
		float P[4][4];			// Estimated Covariance of the state (after measurement z)
		float F[4][4];			// Fundamental Matrix/ State Transition Matrix
		float J[4][6];			// Jacobian of measurement model
		float JRJt[4][4];		// Measurement Quanternion variance
		float SigmaOmega[3]; 	// Gyro spectral noise covariance
		float R[6];				// Measurement noise covariance matrix
		float Q[4][4];
		float dt2;				// half of sampling time
		float magDeclination;	// magnetic declination

	public:
		ExtendedKalmanFilter();
		void SetSampleTime(float freq);
		void SetGyroNoise(float NoiseX, float NoiseY, float NoiseZ);
		void SetR(float NoiseAx, float NoiseAy, float NoiseAz, float NoiseMx, float NoiseMy, float NoiseMz);
		void SetP(float P00, float P11, float P22, float P33);
		bool Run(float ax, float ay, float az, float gx, float gy, float gz, float mx, float my, float mz);
		void GetOrientation(Quaternion& qState);
};

#endif /*EXTENDEDKALMANFILTER_HPP_ */
