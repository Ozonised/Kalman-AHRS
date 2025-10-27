# KALMAN-AHRS
![](https://github.com/Ozonised/Kalman-AHRS/blob/main/testKalman.gif)
An qunaternion based Extended Kalman Filter Attitude Heading and Reference System (AHRS).

## How to use:
1. Include the following source files into your project directory:
  - [qunaterion.hpp](quaternion.hpp)
  - [qunaterion.cpp](quaternion.cpp)
  - [ExtendedKalmanFilter.hpp](ExtendedKalmanFilter.hpp)
  - [ExtendedKalmanFilter.cpp](ExtendedKalmanFilter.cpp)

2. Now create an Quaternion and ExtendedKalmanFilter object.
```cpp
ExtendedKalmanFilter ahrs;
Quaternion q;
```
3. Initialize ```ahrs``` using the following functions: ```ahrs.SetSampleTime(sampling frequency);```
4. Use ```ahrs.Run(accelX, accelY, accelZ, gyroX, gyroY, gyroZ, magX, magY, magZ);``` to compute the orientation quaternion.
5. Use ```ahrs.GetOrientation(q)``` to get the orientation quaternion.

### Additional Functions:
- ```void SetGyroNoise(float NoiseX, float NoiseY, float NoiseZ);``` Sets the gyro noise covariance
- ```void SetR(float NoiseAx, float NoiseAy, float NoiseAz, float NoiseMx, float NoiseMy, float NoiseMz);``` Sets the Measurement Noise Covariance Matrix

## Example:
### Using only Accelerometer and Gyroscope
```cpp
#include "ExtendedKalmanFilter.hpp"
          .
          .
#define IMU_SAMPLING_FREQ 416.0f
ExtendedKalmanFilter ahrs;
Quaternion q;
          .
          .
int main(void)
{
	float accelX, accelY, accelZ, gyroX, gyroY, gyroZ;
	uint8_t ahrsComputeSuccess = 0;
          .
          .
    funcToInitialiseIMU();
    ahrs.SetSampleTime(IMU_SAMPLING_FREQ);

    if (imuDataAvailable()) {
        funcToReadIMUData(accelX, accelY, accelZ, gyroX, gyroY, gyroZ);
        // accelerometer readings should be in m/s^2 and gyro should be in degrees per second
			  ahrsComputeSuccess = ahrs.Run(accelX, accelY, accelZ, gyroX, gyroY, gyroZ, 0, 0,0);
        if (ahrsComputeSuccess) {
            ahrs.GetOrientation(q);
            ahrsComputeSuccess = 0;
        }
    }
}
```
### Note:
When using only accelerometer and gyroscope,the z component of the quaternion is not calculated.

### Using only Accelerometer, Gyroscope and Magnetometer
```cpp
#include "ExtendedKalmanFilter.hpp"
          .
          .
#define IMU_SAMPLING_FREQ 416.0f
ExtendedKalmanFilter ahrs;
Quaternion q;
          .
          .
int main(void)
{
	float accelX, accelY, accelZ, gyroX, gyroY, gyroZ, magX, magY, magZ;
	uint8_t ahrsComputeSuccess = 0;
          .
          .
    funcToInitialiseIMU();
    ahrs.SetSampleTime(IMU_SAMPLING_FREQ);
          .
          .
    if (imuDataAvailable()) {
        funcToReadIMUData(accelX, accelY, accelZ, gyroX, gyroY, gyroZ, magX, magY, magZ);
        // accelerometer readings should be in m/s^2, gyro should be in degrees per second and magnetometer should be in uT
		ahrsComputeSuccess = ahrs.Run(accelX, accelY, accelZ, gyroX, gyroY, gyroZ, magX, magY, magZ);
        if (ahrsComputeSuccess) {
            ahrs.GetOrientation(q);
            ahrsComputeSuccess = 0;
        }
    }
}
```
## Memory Usage
- The ```ahrs``` object consumes around - ```572``` bytes.
- ```bool Run(float ax, float ay, float az, float gx, float gy, float gz, float mx, float my, float mz)``` consumes around - ```352``` bytes.
