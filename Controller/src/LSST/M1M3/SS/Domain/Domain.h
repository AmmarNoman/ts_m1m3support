/*
 * Domain.h
 *
 *  Created on: Nov 2, 2018
 *      Author: ccontaxis
 */

#ifndef LSST_M1M3_SS_DOMAIN_DOMAIN_H_
#define LSST_M1M3_SS_DOMAIN_DOMAIN_H_

int DC_ACCELEROMETER_1X = 0;
int DC_ACCELEROMETER_1Y = 1;
int DC_ACCELEROMETER_2X = 2;
int DC_ACCELEROMETER_2Y = 3;
int DC_ACCELEROMETER_3X = 4;
int DC_ACCELEROMETER_3Y = 5;
int DC_ACCELEROMETER_4X = 6;
int DC_ACCELEROMETER_4Y = 7;
int DC_ACCELEROMETER_ANGULAR_ACCELERATION_X = 0;
int DC_ACCELEROMETER_ANGULAR_ACCELERATION_Y = 1;
int DC_ACCELEROMETER_ANGULAR_ACCELERATION_Z = 2;

struct DCAccelerometerData {
	double Timestamp;
	float RawAccelerometerValues[8];
	float AccelerometerValues[8];
	float AngularAcceleration[3];
};

int IMS_DISPLACEMENT_1 = 0;
int IMS_DISPLACEMENT_2 = 0;
int IMS_DISPLACEMENT_3 = 0;
int IMS_DISPLACEMENT_4 = 0;
int IMS_DISPLACEMENT_5 = 0;
int IMS_DISPLACEMENT_6 = 0;
int IMS_DISPLACEMENT_7 = 0;
int IMS_DISPLACEMENT_8 = 0;
int IMS_POSITION_X = 0;
int IMS_POSITION_Y = 1;
int IMS_POSITION_Z = 2;
int IMS_ROTATION_X = 3;
int IMS_ROTATION_Y = 4;
int IMS_ROTATION_Z = 5;
long IMS_ERROR_UNKNOWN_PROBLEM = 1 << 0;
long IMS_ERROR_INVALID_RESPONSE = 1 << 1;
long IMS_ERROR_RESPONSE_TIMEOUT = 1 << 2;
long IMS_ERROR_SENSOR_REPORTS_INVALID_COMMAND = 1 << 3;
long IMS_ERROR_SENSOR_REPORTS_COMMUNICATION_TIMEOUT_ERROR = 1 << 4;
long IMS_ERROR_SENSOR_REPORTS_DATA_LENGTH_ERROR = 1 << 5;
long IMS_ERROR_SENSOR_REPORTS_NUMBER_OF_PARAMETERS_ERROR = 1 << 6;
long IMS_ERROR_SENSOR_REPORTS_PARAMETER_ERROR = 1 << 7;
long IMS_ERROR_SENSOR_REPORTS_COMMUNICATION_ERROR = 1 << 8;
long IMS_ERROR_SENSOR_REPORTS_ID_NUMBER_ERROR = 1 << 9;
long IMS_ERROR_SENSOR_REPORTS_EXPANSION_LINE_ERROR = 1 << 10;
long IMS_ERROR_SENSOR_REPORTS_WRITE_CONTROL_ERROR = 1 << 11;

struct IMSData {
	double Timestamp;
	float RawDisplacementValues[8];
	float MirrorPosition[6];
	long ErrorFlags;
};


#endif /* LSST_M1M3_SS_DOMAIN_DOMAIN_H_ */
