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
3. Initialize ```ahrs``` using the following functions: ```	ahrs.SetSampleTime(sampling frequency);```
4. Use ```ahrs.Run(accelX, accelY, accelZ, gyroX, gyroY, gyroZ, magX, magY, magZ);``` to compute the orientation quaternion.
5. Use ```ahrs.GetOrientation(q)``` to get the orientation quaternion.
