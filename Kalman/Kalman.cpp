/*
 * Author:        Jacob Morgan
 * Date Created:  07/01/2020
 * Company:       UNCC
 * Extra:         
 * 
 * Description:   This code runs a Kalman filter based on
 *								a sinusoidal model. The user must set what
 *								the model is.
 *
 * How to run:
 * In Setup:
 *	1. initialize()
 *	2. setWheelCircumference_m(float wheelCircumference_m_);
 * 	3. setStdDevModel(float stdDev);
 *	4. setStdDevSensor(float stdDev);
 *	5. setSinusoidBounds(float bound);
 *
 * When turning on motor to turn wheel
 *	6. setWheelPeriod_s(float period_s);
 *  7. setWheelDirection(direction wheelDirection_);
 *  8. setMaxMinSensorValues(float min_, float max_);
 *	9. initKalmanPrediction(float predictionX, float predictionY);
 *
 *	Then can turn on motors after doing above steps and entering the loop
 *
 * Then in loop:
 *	10.setSensorReadTimeDelta(float timeDelta_s);
 *	11.setSensorReadings(float sensorReading_mpss);
 *	12.runFilter();
 *	13.get Kalman Results using getPredictionValueX() and getPredictionValueY()
 *  14.get distance traveled using getDistanceTaveled_m()
 *
*/

#include <Kalman.h>
#include <math.h>
#include <HardwareSerial.h>
#include <WString.h>

#define constrain(amt,low,high) ((amt)<(low)?(low):((amt)>(high)?(high):(amt)))
#define SQUARESIZE 3

// Constructor
Kalman::Kalman() {}
// De-Constructor
Kalman::~Kalman(){}

void Kalman::initialize() {
	uint8_t row, col;
	for(row = 0; row < SQUARESIZE; row++) {
		for(col = 0; col < SQUARESIZE; col++) {
			if(row == 1 && col == 1) {								// Initialize derivative of the state transition equations and covariance of the measurement noises
				dfda[row][col] = 1.0;
				Q[row][col] = 1.0;
			} else {
				dfda[row][col] = 0.0;
				Q[row][col] = 0.0;
			}
			
			if(row == 0 && col == 0) {								// Initialize derivative of the state transition equations
				dfdx[row][col] = 1.0;
			} else if(row == 0 && col == 1) {
				dfdx[row][col] = 1.0;
			} else if(row == 1 && col == 1) {
				dfdx[row][col] = 1.0;
			} else if(row == 2 && col == 0) {
				dfdx[row][col] = 1.0;
			} else {
				dfdx[row][col] = 0.0;
			}
			
			if(row == col) {													// Initialize identity matrix and Covariance
				eye[row][col] = 1.0;
				S[row][col] = 1.0;
				R[row][col] = 1.0;
			} else {
				eye[row][col] = 0.0;
				S[row][col] = 0.0;
				R[row][col] = 0.0;
			}
		}	// end col for loop
	}	// end row for loop
	
	X[0] = 0.0;
	X[1] = 0.0;
	X[2] = 0.0;
	
	Matrix.Print((mtx_type*)eye, 3, 3, "I"); 
}

void Kalman::runFilter() {
	if(abort) {
		return;
	}
	
	if(!checkSetup) {
		if(!checkSetupComplete()) {
			Serial.println(F("\nDid not complete Setup. Ensure that all values are set"));
			while(1){}
		}
	}
	
	//----------------------------------------------------------
	//-----------Predict----------------------------------------
	//----------------------------------------------------------
	
	//------------------------------
	// prediction for state X
	//------------------------------
	lastPrediction = X[2];
	if(firstRunDone) {
		X[0] = fmod((X[0] + X[1]*deltaT), 2*PI);	// Keep between 0 and 2*PI
		// X[1] will not change in the prediction step. It is expected to stay constant.
	} else {
		/* Serial.println(F("In First Iteration")); */
		X[0] = fmod((2*PI*(1/wheelPeriod_s)*modelTime_s), 2*PI);
		X[1] = 2*PI*(1/wheelPeriod_s);
		Serial.print(F("Initial dx/dt = ")); Serial.println(X[1], 6);
	}
	
	X[2] = amplitude*sin(X[0]) - ampShift;
	
	//------------------------------
	// prediction for Covariance		// Only update covariance, S, after first kalman iteration
	//------------------------------
	if(firstRunDone) {
		dfdx[0][1] = deltaT;
		dfdx[2][0] = amplitude*cos(X[0]);
		
		mtx_type dfdxT[3][3];
		Matrix.Transpose((mtx_type*) dfdx, 3, 3, (mtx_type*) dfdxT);	// df/dx'
		
		mtx_type dfdx_S[3][3];
		Matrix.Multiply((mtx_type*) dfdx, (mtx_type*) S, 3, 3, 3, (mtx_type*) dfdx_S);	// dfdx*S
		
		mtx_type dfdx_S_dfdxT[3][3];
		Matrix.Multiply((mtx_type*) dfdx_S, (mtx_type*) dfdxT, 3, 3, 3, (mtx_type*) dfdx_S_dfdxT);	// dfdx*S*dfdx'
		
		// dfda transpose is just dfda since a scalar value sits on the diagonal and zero everywhere else
		// Since dfda has 1 in center and zero everywhere else, dfda*Q*dfda' = Q
		Matrix.Add((mtx_type*) dfdx_S_dfdxT, (mtx_type*) Q, 3, 3, (mtx_type*) S);	// S = dfdx*S*dfdx' + dfda*Q*dfda'
	}
	
	//----------------------------------------------------------
	//-----------Update-----------------------------------------
	//----------------------------------------------------------
	
	//------------------------------
	// Determine K, kalman gain
	//------------------------------
	mtx_type SR[3][3];
	Matrix.Add((mtx_type*) S, (mtx_type*) R, 3, 3, (mtx_type*) SR);	// S + R
	
	if(!Matrix.Invert((mtx_type*) SR, 3)) {		// (S + R)^-1
		Serial.println(F("Failed to Invert Matrix. Skipping Update"));
		return;
	}
	
	mtx_type K[3][3];	// Kalman Gain
	Matrix.Multiply((mtx_type*) S, (mtx_type*) SR, 3, 3, 3, (mtx_type*) K);	// S * (S + R)^-1
	
	//------------------------------
	// Determine updated X
	//------------------------------
	mtx_type YmX[3];
	Matrix.Subtract((mtx_type*) sensorModel, (mtx_type*) X, 3, 1, (mtx_type*) YmX);				//        (Y - X)
	
	mtx_type K_YmX[3];
	Matrix.Multiply((mtx_type*) K, (mtx_type*) YmX, 3, 3, 1, (mtx_type*) K_YmX);					//     K *(Y - X)
	
	mtx_type newX[3];
	Matrix.Add((mtx_type*) X, (mtx_type*) K_YmX, 3, 1, (mtx_type*) newX);									// X + K *(Y - X), 3x1
	
	X[0] = newX[0];
	X[1] = newX[1];
	X[2] = newX[2];
	
	//------------------------------
	// Determine updated Covariance			// K*dgdx = K
	//------------------------------
	mtx_type ImK_dgdx[3][3];
	Matrix.Subtract((mtx_type*) eye, (mtx_type*) K, 3, 3, (mtx_type*) ImK_dgdx);			//  I - K*dgdx
	
	mtx_type newS[3][3];
	Matrix.Multiply((mtx_type*) ImK_dgdx, (mtx_type*) S, 3, 3, 3, (mtx_type*) newS);			// [I - K*dgdx]*S
	
	uint8_t row, col;
	for(row = 0; row < SQUARESIZE; row++) {
		for(col = 0; col < SQUARESIZE; col++) {
			S[row][col] = newS[row][col];
		}
	}
	
	//----------------------------------------------------------
	//-----------Get Distance-----------------------------------
	//----------------------------------------------------------
	
	// Determine distance traveled with new results
	if(firstRunDone) {
		getDistanceTraveled();
	} else {
		firstRunDone = true;
	}
}

void Kalman::initKalmanPrediction(float predictionX, float predictionY) {
	if(!maxMinSensorValuesSet) {
		Serial.println(F("Aborting application. Need to set wheel period, circumference, and direction before init kalman filter"));
		Serial.println(F("Need max and min sensor values set"));
		while(1){}
	}
	
	X[2] = predictionX;
	lastPrediction = predictionX;
	
	float velocity = wheelCircumference_m/wheelPeriod_s;
	float phi	= (predictionX + ampShift)/calculatedAmplitude;
	if(phi > 1) {
		phi = 1;
	} else if(phi < -1) {
		phi = -1;
	}
	
	modelTime_s = asin(phi)*wheelPeriod_s/(2*PI);
	
	// ASIN will only return a value in the 1st or 4th quadrantt. These Statement appropriately
	//	account for time offset if sine should be in the 2nd or 3rd quadrant. Starting Quadrant
	//	is determined by sign of starting x and starting y.
	if(predictionX >= 0) {
		if(predictionY >= 0) {		// X is in quadrant 1, y in 2
			currentQuadrant = 1;
		} else {									// X is in quadrant 2, y in 3
			currentQuadrant = 2;
			modelTime_s = wheelPeriod_s/2 - modelTime_s;
		}
	} else {
		if(predictionY >= 0) {		// X is in quadrant 4, y in 1
			currentQuadrant = 4;
			modelTime_s += wheelPeriod_s;
		} else {									// X is in quadrant 3, y in 4
			currentQuadrant = 3;
			modelTime_s = fabs(modelTime_s) + wheelPeriod_s/2;
		}
	}
	
	kalmanIsInit = true;
}

void Kalman::setWheelDirection(direction wheelDirection_) {
	wheelDirection = wheelDirection_;
	wheelDirectionSet = true;
}

void Kalman::setWheelCircumference_m(float wheelCircumference_m_) {
	wheelCircumference_m = wheelCircumference_m_;
	wheelCircumferenceSet = true;
	
	/* Serial.print(F("Wheel Circumference = ")); Serial.println(wheelCircumference_m_); */
}

void Kalman::setWheelPeriod_s(float period_s) {
	wheelPeriod_s 	= period_s;
	wheelPeriodSet 	= true;
}

void Kalman::setStdDevModel(float stdDev) {
	std_dev_model		= stdDev;
	Q[1][1] = pow(stdDev, 2);
	stdDevModelSet 	= true;
}

void Kalman::setStdDevSensor(float stdDev) {
	std_dev_sensor	= stdDev;
	R[0][0] = pow(std_dev_sensor, 2);
	R[1][1] = pow(std_dev_sensor, 2);
	R[2][2] = pow(std_dev_sensor, 2);
	stdDevSensorSet = true;
}

void Kalman::setMaxMinSensorValues(float min_, float max_) {
	if(!wheelPeriodSet || !wheelDirectionSet || !wheelCircumferenceSet)
	{
		Serial.println(F("Aborting application. Need to set wheel period, circumference, and direction before setting max/min"));
		while(1){}
	}
	
	ampShift = (fabs(min_) - fabs(max_))/2;
	calculatedAmplitude = ((max_ + ampShift) - (-min_ + ampShift))/2;
	float velocity = wheelCircumference_m/wheelPeriod_s;
	sensorCenterOffset = pow(velocity, 2)/(calculatedAmplitude - GRAVITYMAGNITUDE);
	centrifugal = pow(velocity, 2)/sensorCenterOffset;
	amplitude = GRAVITYMAGNITUDE + centrifugal;
	maxMinSensorValuesSet	= true;
	
	/* Serial.print(F("Offset = ")); Serial.println(ampShift, 6);
	Serial.print(F("Amplitude = ")); Serial.println(calculatedAmplitude, 6);
	Serial.print(F("Velocity = ")); Serial.println(velocity, 6);
	Serial.print(F("Sensor Center Offset = ")); Serial.println(sensorCenterOffset, 6);
	Serial.print(F("centrifugal = ")); Serial.println(centrifugal, 6); */
}

void Kalman::setSinusoidBounds(float bound) {
	bounds = bound;
	boundsSet = true;
}

#define AVERAGE 5
void Kalman::setSensorReadings(float sensorReading_mpss) {
	static bool averageBufferFull = false;
	static uint8_t index = 0;
	static float averageBuffer[AVERAGE] = {0};
	
	averageBuffer[index] = sensorReading_mpss;
	
	if(!averageBufferFull && index == AVERAGE - 1) {
		averageBufferFull = true;
	}
	
	index++;
	index %= AVERAGE;
	
	uint8_t i = 0;
	float sum = 0;
	for(i = 0; i < AVERAGE; i++) {
		sum += averageBuffer[i];
		//Serial.print(F("Index ")); Serial.print(i); Serial.print(F(" = ")); Serial.println(averageBuffer[i], 6);
	}
	
	if(averageBufferFull) {
		sum /= AVERAGE;
	} else {
		sum /= index;
	}
	//Serial.print(F("Sum = ")); Serial.println(sum, 6);
	
	if(firstRunDone) {
		sensorModel[0] = fmod((sensorModel[0] + sensorModel[1]*deltaT), 2*PI);	// Keep between 0 and 2*PI
	} else {
		sensorModel[0] = fmod((2*PI*(1/wheelPeriod_s)*modelTime_s), 2*PI);
	}
	sensorModel[1] = 2 * PI * (1 / wheelPeriod_s);
	sensorModel[2] = sum;
	
	/* Serial.print(F("Sensor Angle = ")); Serial.println(sensorModel[0], 6);
	Serial.print(F("Sensor dx/dt = ")); Serial.println(sensorModel[1], 6);
	Serial.print(F("Sensor accel = ")); Serial.println(sensorModel[2], 6); */
}

void Kalman::setSensorReadTimeDelta(float timeDelta_s) {
	deltaT = timeDelta_s;
}

bool Kalman::checkSetupComplete() {
	checkSetup = false;
	if(kalmanIsInit	&& wheelPeriodSet && stdDevModelSet && stdDevSensorSet && boundsSet && wheelDirectionSet && wheelCircumferenceSet && maxMinSensorValuesSet) {
		checkSetup = true;
	}
	return checkSetup;
}

void Kalman::resetKalmanFilter() {
	checkSetup 						= false;		// Verify Set up was checked 
	kalmanIsInit					= false;		// Checks model time is set
	wheelPeriodSet				= false;
	wheelDirectionSet			= false;
	maxMinSensorValuesSet	= false;
	modelTime_s						= 0;
	X[2] 									= 0;
	// Will Keep stdDevModel the same
	// Will Keep stdDevSensor the same
	// Will Keep stdDevDynamic the same
	// Will Keep bounds the same
}

void Kalman::getDistanceTraveled() {
	float radToDistance_m = wheelCircumference_m/(2*PI);
	float phi	= (X[2] + ampShift)/amplitude;
	if(phi > 1) {
		phi = 1;
	} else if(phi < -1) {
		phi = -1;
	}
	
	calculatedAngle_rad = asin(phi);
	
	// Determine what quadrant sinusoid is in
	uint8_t lastQuadrant = currentQuadrant;
	currentQuadrant = getNextQuadrant(X[2]);
	calculatedAngle_rad = getCalculatedAngle_rad(currentQuadrant, calculatedAngle_rad);
	
	
	if(lastQuadrant == 4 && currentQuadrant == 1) {
		distanceTraveled_m += calculatedAngle_rad*radToDistance_m;
	} else {
		distanceTraveled_m += (calculatedAngle_rad - lastCalculatedAngle_rad)*radToDistance_m;
	}
	lastCalculatedAngle_rad = calculatedAngle_rad;
}

uint8_t Kalman::getNextQuadrant(float currentAcceleration_mpss) {
	uint8_t quadrant;
	if(currentQuadrant == 1) {
		if(currentAcceleration_mpss < lastPrediction) {
			quadrant = 2;
		}
	} else if(currentQuadrant == 2) {
		if(currentAcceleration_mpss < 0) {
			quadrant = 3;
		}
	} else if(currentQuadrant == 3) {
		if(currentAcceleration_mpss > lastPrediction) {
			quadrant = 4;
		}
	} else {
		if(currentAcceleration_mpss > 0) {
			quadrant = 1;
		}
	}
	
	return quadrant;
}

float Kalman::getCalculatedAngle_rad(uint8_t quadrant, float currentAngle_rad) {
	// Determine calculated angle in radians
	if(quadrant == 1) {
		currentAngle_rad = fabs(currentAngle_rad);
	} else if(quadrant == 2) {
		currentAngle_rad = PI - currentAngle_rad;
	} else if(quadrant == 3) {
		currentAngle_rad = fabs(currentAngle_rad) + PI;
	} else {
		currentAngle_rad += 2*PI;
	}
	
	return currentAngle_rad;
}

//	End of File
