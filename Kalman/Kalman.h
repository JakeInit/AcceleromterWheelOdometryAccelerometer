/*
 * Author:        Jacob Morgan
 * Date Created:  07/01/2020
 * Company:       UNCC
 * Extra:         
 * 
 * Description:   This code runs a Kalman filter based on
 *								a sinusoidal model. The user must set what
 *								the model is.
*/

#ifndef _KALMAN_H_
#define _KALMAN_H_

#include <stdint.h>

#ifndef PI
#define PI  3.1415926535897932384626433832795
#endif

#define GRAVITYMAGNITUDE	9.80665

template <typename T> int sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}

class Kalman
{
public:
	Kalman();
	~Kalman();
	
	enum direction{CLOCKWISE, COUNTERCLOCK};

	void 	runFilter();
	void	initialize();
	void 	initKalmanPrediction(float predictionX, float predictionY);
	void	setWheelDirection(direction wheelDirection_);
	void	setWheelCircumference_m(float wheelCircumference_m_);
	void 	setWheelPeriod_s(float period_s);
	void 	setStdDevModel(float stdDev);
	void 	setStdDevSensor(float stdDev);
	void	setMaxMinSensorValues(float min_, float max_);
	void	setSinusoidBounds(float bound);
	void	setSensorReadings(float sensorReading_mpss);
	void	setSensorReadTimeDelta(float timeDelta_s);
	void	resetKalmanFilter();
	
private:
	void	determineDistanceTraveled();
	bool	checkSetupComplete();
	
	bool abort 									= false;
	bool kalmanIsInit						= false;
	bool firstRunDone						= false;
	bool wheelPeriodSet					= false;
	bool stdDevModelSet					= false;
	bool stdDevSensorSet				= false;
	bool maxMinSensorValuesSet	= false;
	bool boundsSet							= false;
	bool wheelDirectionSet			= false;
	bool wheelCircumferenceSet	= false;
	bool checkSetup							= false;
	
	direction wheelDirection = CLOCKWISE;
	
	uint8_t quadrant	= 1;		// Quadrant of unit circle that the sinusoid is in
	
	float modelTime_s = 0;
	float sensorValue = 0;
	float lastPrediction = 0;
	float calculatedAngle_rad = 0;
	float lastCalculatedAngle_rad = 0;
	float amplitude = 0;
	
	float centrifugal	= 0;
	float ampShift 		= 0;
	float calculatedAmplitude = 0;	// Amplitude of centered sinewave. Sinewave does not center on zero since sensor not perfectly on middle of wheel
	float sensorCenterOffset 	= 0;	// Sensor Offset from center of wheel
	
	// Caclculated Every Interval and will have initial values in Constructor
	float X[3];			// State transition matrix
	float S[9];			// Covariance Matrix
	float dfdx[9];	// Derivative of the state transition equations with respect to the state variables
									//		T (position 1) gets set in setSensorReadTimeDelta();
	
	// Need initialized by user and will have initial value in Constructor
	float Q[9];			// Covariance of the dynamic noises				is set in setStdDevModel
	float R;				// Covariance of the measurement noises		is set in setStdDevSensor
	
	// Needs initialized just in Constructor
	float eye[9];		// 3x3 identity matrix
	float dfda[9];	// Derivative of the state transition equations with respect to the dynamic noises
	float dgdx[3];	// Derivative of the observation equations with respect to the state variables
	
	// initialized as constant
	// const float dgdn = 1;	// Derivative of the observation equations with respect to the measurement noises
	
	
	float wheelPeriod_s		= 0;
	float std_dev_model		= 0;
	float std_dev_sensor	= 0;
	float std_dev_dynamic	= 0;
	float bounds					= 1;
	float deltaT 					= 0;				// Time between measurements
	float wheelCircumference_m	= 0;
	double distanceTraveled_m		= 0;	// Can be based on x or y axis or by using both
	
public:
	bool	getAbortStatus()				const {return abort;}
	bool	getKalmanInitStatus()		const {return kalmanIsInit;}
	float getPredictionValueX()		const {return X[2];}
	float getDistanceTaveled_m()	const {return distanceTraveled_m;}
};

#endif