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

#define FLOAT sizeof(float)
//#define DEBUG

void multiplyMatrixAxB(const uint8_t rows1, const uint8_t cols1, float size1, float matrix1[], const uint8_t rows2, const uint8_t cols2, float size2, float matrix2[]);
void addMatrixAnB(const uint8_t rows1, const uint8_t cols1, float size1, float matrix1[], const uint8_t rows2, const uint8_t cols2, float size2, float matrix2[]);
void subMatrixAnB(const uint8_t rows1, const uint8_t cols1, float size1, float matrix1[], const uint8_t rows2, const uint8_t cols2, float size2, float matrix2[]);
void multiplyMatrixByScalar(const uint8_t rows_, const uint8_t cols_, float size_, float matrix[], const float scalar);
void transposeMatrix(const uint8_t rows_, const uint8_t cols_, float size_, float matrix[]);

void debugPrompt();
void print3x3(float matrix[]);
void print3x1(float matrix[]);
bool checkNaN3x1(float matrix[]);
bool checkNaN3x3(float matrix[]);
// void inverseMatrix(const uint8_t order, float matrix[]);

// Temp Variable to hold matrixes
float temp3x3[9] = {0};
float temp3[3] = {0};

// Determines if size of array is 9 or 3
bool squared = false;

// Constructor
Kalman::Kalman() {}
// De-Constructor
Kalman::~Kalman(){}

void Kalman::initialize() {
	uint8_t i;
	for(i = 0; i < 9; i++)
	{
		if(i == 4)																	// Initialize derivative of the state transition equations and covariance of the measurement noises
		{
			dfda[i] = 1.0;
			Q[i] = 1.0;
		}
		else
		{
			dfda[i] = 0.0;
			Q[i] = 0.0;
		}
		
		if(i == 0 || i == 1 || i == 4 || i == 6)		// Initialize derivative of the state transition equations
		{
			dfdx[i] = 1.0;
		}
		else
		{
			dfdx[i] = 0.0;
		}
		
		if(i == 0 || i == 4 || i == 8)							// Initialize identity matrix and Covariance
		{
			eye[i] = 1.0;
			S[i] = 1.0;
			R[i] = 1.0;
		}
		else
		{
			eye[i] = 0.0;
			S[i] = 0.0;
			R[i] = 0.0;
		}
	}
	
	X[0] = 0.0;
	X[1] = 0.0;
	X[2] = 0.0;
	
	/* print3x1(X);
	print3x3(dgdx);
	print3x3(dfda);
	print3x3(dfdx);
	print3x3(S); */
}

void Kalman::runFilter() {
	if(abort)
	{
		return;
	}
	
	if(!checkSetup)
	{
		if(!checkSetupComplete())
		{
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
	if(firstRunDone)
	{
		/* Serial.print(F("X[0] = ")); Serial.print(X[0], 6); */
		X[0] = fmod((X[0] + X[1]*deltaT), 2*PI);	// Keep between 0 and 2*PI
		// X[1] will not change in the prediction step. It is expected to stay constant.
	}
	else
	{
		/* Serial.println(F("In First Iteration")); */
		X[0] = fmod((2*PI*(1/wheelPeriod_s)*modelTime_s), 2*PI);
		X[1] = 2*PI*(1/wheelPeriod_s);
	}
	
	X[2] = amplitude*sin(X[0]) - ampShift;
	
	#ifdef DEBUG
		Serial.print(F("X = "));
		print3x1(X);
		debugPrompt();
	#endif
	
	//------------------------------
	// prediction for Covariance
	//------------------------------
	if(firstRunDone)										// Only update covariance, S, after first kalman iteration
	{
		dfdx[1] = deltaT;
		dfdx[6] = amplitude*cos(X[0]);
		
		squared = true;
		float dfdxTranspose[9];
		transposeMatrix(3, 3, sizeof(dfdx)/FLOAT, dfdx);
		memcpy(dfdxTranspose, temp3x3, sizeof(temp3x3));		// df/dx'
		
		multiplyMatrixAxB(3, 3, sizeof(dfdx)/FLOAT, dfdx, 3, 3, sizeof(S)/FLOAT, S);																// dfdx*S
		multiplyMatrixAxB(3, 3, sizeof(temp3x3)/FLOAT, temp3x3, 3, 3, sizeof(dfdxTranspose)/FLOAT, dfdxTranspose);	// dfdx*S*dfdx'
		
		/* if(checkNaN3x3(temp3x3))
		{
			Serial.println(F("Abort Program due to nan"));
			Serial.print(F("left Covariance ="));
			print3x3(temp3x3);
			Serial.print(F("dfdx ="));
			print3x3(dfdx);
			Serial.print(F("dfdxTranspose ="));
			print3x3(dfdxTranspose);
			Serial.print(F("S = "));
			print3x3(S);
			abort = true;
			return;
		} */
		
		float covarianceLeft[9];
		memcpy(covarianceLeft, temp3x3, sizeof(temp3x3));
		
		// dfda transpose is just dfda since a scalar value sits on the diagonal and zero everywhere else
		// Since dfda has 1 in center and zero everywhere else, dfda*Q*dfda' = Q
		
		addMatrixAnB(3, 3, sizeof(covarianceLeft)/FLOAT, covarianceLeft, 3, 3, sizeof(Q)/FLOAT, Q);	// dfdx*S*dfdx' + dfda*Q*dfda'
		
		memcpy(S, temp3x3, sizeof(temp3x3));	// predicted covariance, S, 3x3
	}
	
	#ifdef DEBUG
		Serial.print(F("S = "));
		print3x3(S);
		debugPrompt();
	#endif
	
	//----------------------------------------------------------
	//-----------Update-----------------------------------------
	//----------------------------------------------------------
	
	//------------------------------
	// Determine K, kalman gain
	//------------------------------
	squared = true;
	addMatrixAnB(3, 3, sizeof(S)/FLOAT, S, 3, 3, sizeof(R)/FLOAT, R);	// S + R
	multiplyMatrixAxB(3,3, sizeof(S)/FLOAT, S, 3, 3, sizeof(temp3x3)/FLOAT, temp3x3);	// S * (S + R)
	
	float K[9];														// This is the Kalman Gain (3x3)
	memcpy(K, temp3x3, sizeof(temp3x3));
	
	#ifdef DEBUG
		Serial.print(F("K = "));
		print3x3(K);
		debugPrompt();
	#endif
	
	//------------------------------
	// Determine updated X
	//------------------------------
	squared = false;
	subMatrixAnB(3, 1, sizeof(sensorModel)/FLOAT, sensorModel, 3, 1, sizeof(X)/FLOAT, X);	//        (Y - X)
	multiplyMatrixAxB(3, 3, sizeof(K)/FLOAT, K, 3, 1, sizeof(temp3)/FLOAT, temp3);				//     K *(Y - X)
	addMatrixAnB(3, 1, sizeof(X)/FLOAT, X, 3, 1, sizeof(temp3)/FLOAT, temp3);							// X + K *(Y - X), 3x1
	memcpy(X, temp3, sizeof(temp3));																											// This updates state matrix, X, 3x1
	
	#ifdef DEBUG
		Serial.print(F("Updated X = "));
		print3x1(X);
		debugPrompt();
	#endif
	
	//------------------------------
	// Determine updated Covariance
	//------------------------------
	/* Serial.print(F("K")); print3x1(K); */
	// K*dgdx = K
	squared = true;
	subMatrixAnB(3, 3, sizeof(eye)/FLOAT, eye, 3, 3, sizeof(K)/FLOAT, K);							//  I - K*dgdx
	multiplyMatrixAxB(3,3, sizeof(temp3x3)/FLOAT, temp3x3, 3, 3, sizeof(S)/FLOAT, S);	// [I - K*dgdx]*S
	memcpy(S, temp3x3, sizeof(temp3x3));	// This is updated covariance, S, 3x3
	
	#ifdef DEBUG
		Serial.print(F("Updated S = "));
		print3x3(S);
		debugPrompt();
	#endif
	
	//----------------------------------------------------------
	//-----------Get Distance-----------------------------------
	//----------------------------------------------------------
	
	// Determine distance traveled with new results
	if(firstRunDone)
	{
		getDistanceTraveled();
	}
	else
	{
		firstRunDone = true;
	}
}

void Kalman::initKalmanPrediction(float predictionX, float predictionY) {
	if(!maxMinSensorValuesSet)
	{
		Serial.println(F("Aborting application. Need to set wheel period, circumference, and direction before init kalman filter"));
		Serial.println(F("Need max and min sensor values set"));
		while(1){}
	}
	
	X[2] = predictionX;
	lastPrediction = predictionX;
	
	float velocity = wheelCircumference_m/wheelPeriod_s;
	float phi	= (predictionX + ampShift)/calculatedAmplitude;
	if(phi > 1)
	{
		phi = 1;
	}
	else if(phi < -1)
	{
		phi = -1;
	}
	
	modelTime_s = asin(phi)*wheelPeriod_s/(2*PI);
	
	// ASIN will only return a value in the 1st or 4th quadrantt. These Statement appropriately
	//	account for time offset if sine should be in the 2nd or 3rd quadrant. Starting Quadrant
	//	is determined by sign of starting x and starting y.
	if(predictionX >= 0)
	{
		if(predictionY >= 0)			// X is in quadrant 1, y in 2
		{
			currentQuadrant = 1;
		}
		else											// X is in quadrant 2, y in 3
		{
			currentQuadrant = 2;
			modelTime_s = wheelPeriod_s/2 - modelTime_s;
		}
	}
	else
	{
		if(predictionY >= 0)			// X is in quadrant 4, y in 1
		{
			currentQuadrant = 4;
			modelTime_s += wheelPeriod_s;
		}
		else											// X is in quadrant 3, y in 4
		{
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
	Q[4] = pow(stdDev, 2);
	stdDevModelSet 	= true;
	
	/* Serial.print(F("Q = "));
	print3x3(Q); */
}

void Kalman::setStdDevSensor(float stdDev) {
	std_dev_sensor	= stdDev;
	R[0] = pow(std_dev_sensor, 2);
	R[4] = pow(std_dev_sensor, 2);
	R[8] = pow(std_dev_sensor, 2);
	stdDevSensorSet = true;
	
	/* Serial.print(F("R = ")); Serial.println(R, 6); */
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

void Kalman::setSensorReadings(float sensorReading_mpss) {
	float phi = (sensorReading_mpss - ampShift)/amplitude;
	if(phi > 1) {
		phi = 1;
	} else if(phi < -1) {
		phi = -1;
	}
	
	float sensorAngle = asin(phi);
	uint8_t quadrant = getNextQuadrant(sensorReading_mpss);
	sensorAngle = getCalculatedAngle_rad(quadrant, sensorAngle);
	
	sensorModel[0] = sensorAngle;
	sensorModel[1] = (sensorAngle - lastCalculatedAngle_rad)*(1/deltaT);
	sensorModel[2] = sensorReading_mpss;
}

void Kalman::setSensorReadTimeDelta(float timeDelta_s) {
	deltaT = timeDelta_s;
}

bool Kalman::checkSetupComplete() {
	checkSetup = false;
	if(kalmanIsInit	&& wheelPeriodSet && stdDevModelSet && stdDevSensorSet && boundsSet && wheelDirectionSet && wheelCircumferenceSet && maxMinSensorValuesSet)
	{
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
	
	lastCalculatedAngle_rad = calculatedAngle_rad;
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
	
}

uint8_t Kalman::getNextQuadrant(float currentAcceleration_mpss) {
	uint8_t quadrant;
	if(currentQuadrant == 1)
	{
		if(currentAcceleration_mpss < lastPrediction)
		{
			quadrant = 2;
		}
	}
	else if(currentQuadrant == 2)
	{
		if(currentAcceleration_mpss < 0)
		{
			quadrant = 3;
		}
	}
	else if(currentQuadrant == 3)
	{
		if(currentAcceleration_mpss > lastPrediction)
		{
			quadrant = 4;
		}
	}
	else
	{
		if(currentAcceleration_mpss > 0)
		{
			quadrant = 1;
		}
	}
	
	return quadrant;
}

float Kalman::getCalculatedAngle_rad(uint8_t quadrant, float currentAngle_rad) {
	// Determine calculated angle in radians
	if(quadrant == 1)
	{
		currentAngle_rad = fabs(currentAngle_rad);
	}
	else if(quadrant == 2)
	{
		currentAngle_rad = PI - currentAngle_rad;
	}
	else if(quadrant == 3)
	{
		currentAngle_rad = fabs(currentAngle_rad) + PI;
	}
	else
	{
		currentAngle_rad += 2*PI;
	}
	
	return currentAngle_rad;
}





//---------------------------------------------------------------------
//---------Static Functions--------------------------------------------
//---------------------------------------------------------------------

void multiplyMatrixAxB(const uint8_t rows1, const uint8_t cols1, float size1, float matrix1[], const uint8_t rows2, const uint8_t cols2, float size2, float matrix2[]) {
  if(cols1 != rows2)
  {
    Serial.println(F("Multiplying Invalid Matrixes. Check ColsA and RowsB"));
    while(1){}
  }
  
  if(size1 != rows1*cols1)
  {
    Serial.println(F("rows and columns do not match size of matrix A"));
    while(1){}
  }
  
  if(size2 != rows2*cols2)
  {
    Serial.println(F("rows and columns do not match size of matrix B"));
    while(1){}
  }
  
  uint8_t row, column, k, i;
  float mat1[rows1][cols1] = {{0},{0}};
  float mat2[rows2][cols2] = {{0},{0}};
  float newMat[rows1][cols2] = {{0},{0}};

  // place matrix1 into 2d array
  row = 0;
  for(row = 0; row < rows1; row++)
  {
    column = 0;
    for(column = 0; column < cols1; column++)
    {
			mat1[row][column] = 0.0;
      mat1[row][column] = matrix1[row*cols1 + column];
    }
  }

  // place matrix2 into 2d array
  row = 0;
  for(row = 0; row < rows2; row++)
  {
    column = 0;
    for(column = 0; column < cols2; column++)
    {
			mat2[row][column] = 0.0;
      mat2[row][column] = matrix2[row*cols2 + column];
    }
  }

  // get new array with multiplied values
  row = 0;
  for(row = 0; row < rows1; row++)
  {
    column = 0;
    for(column = 0; column < cols2; column++)
    {
      k = 0;
			newMat[row][column] = 0.0;
      for(k = 0; k < cols1; k++)
      {
        newMat[row][column] += (float) (mat1[row][k]) * (float) (mat2[k][column]);
      }
    }
  }

  row = 0;
  i = 0;
  // Place new matrix into single array format
  for(row = 0; row < rows1; row++)
  {
    column = 0;
    for(column = 0; column < cols2; column++)
    {
      if(squared)
      {
				temp3x3[i] = 0.0;
        temp3x3[i] = newMat[row][column];
      }
      else
      {
				temp3[i] = 0.0;
        temp3[i] = newMat[row][column];
      }
      i++;
    }
  }
}

void addMatrixAnB(const uint8_t rows1, const uint8_t cols1, float size1, float matrix1[], const uint8_t rows2, const uint8_t cols2, float size2, float matrix2[]) {
  if(rows1 != rows2 || cols1 != cols2)
  {
    Serial.println(F("Matrix Addition Failed: rows and columns sizes must match for each matrix"));
    while(1){}
  }
  
  uint8_t i;
  for(i = 0; i < size1; i++)
  {
		if(squared)
		{
			temp3x3[i] = matrix1[i] + matrix2[i];
		}
		else
		{
			temp3[i] = matrix1[i] + matrix2[i];
		}
  }
}

void subMatrixAnB(const uint8_t rows1, const uint8_t cols1, float size1, float matrix1[], const uint8_t rows2, const uint8_t cols2, float size2, float matrix2[]) {
	if(rows1 != rows2 || cols1 != cols2)
  {
    Serial.println(F("Matrix Subtraction Failed: rows and columns sizes must match for each matrix"));
    while(1){}
  }
  
  uint8_t i;
  for(i = 0; i < size1; i++)
  {
		if(squared)
		{
			temp3x3[i] = matrix1[i] - matrix2[i];
		}
		else
		{
			temp3[i] = matrix1[i] - matrix2[i];
		}
  }
}

void multiplyMatrixByScalar(const uint8_t rows_, const uint8_t cols_, float size_, float matrix[], const float scalar) {
  if(size_ != rows_*cols_)
  {
    Serial.println(F("Scalar Multiply Failed: rows and columns don't match matrix size"));
    while(1){}
  }
  
  uint8_t i;
  for(i = 0; i < size_; i++)
  {
		if(squared)
		{
			temp3x3[i] = matrix[i]*scalar;
		}
		else
		{
			temp3[i] = matrix[i]*scalar;
		}
  }
}

void transposeMatrix(const uint8_t rows_, const uint8_t cols_, float size_, float matrix[]) {
  if(size_ != rows_*cols_)
  {
    Serial.println(F("Fail transpose because rows and columns do not match size of matrix"));
    while(1){}
  }

  float mat[rows_][cols_] = {{0}, {0}};
  float newMat[rows_][cols_] = {{0}, {0}};
  uint8_t row, column, i;
  
  // Place matrix of type vector into an array
  row = 0;
  for(row = 0; row < rows_; row++)
  {
    column = 0;
    for(column = 0; column < cols_; column++)
    {
      mat[row][column] = matrix[row*rows_ + column];
    }
  }

  // Implement Transpose
  column = 0;
  for(column = 0; column < cols_; column++)
  {
    row = 0;
    for(row = 0; row < rows_; row++)
    {
      newMat[row][column] = mat[column][row];
    }
  }
  
  row = 0;
  i = 0;
  // Place new matrix into single array format
  for(row = 0; row < rows_; row++)
  {
    column = 0;
    for(column = 0; column < cols_; column++)
    {
      if(squared)
      {
        temp3x3[i] = newMat[row][column];
      }
      else
      {
        temp3[i] = newMat[row][column];
      }
      i++;
    }
  }
}

/* void inverseMatrix3x3(const uint8_t order, float matrix[])
{
  if(9 != pow(order,2))
  {
    Serial.println("Fail Inverse because rows and columns do not match size of matrix");
    while(1){}
  }
	
  float determinant = 0;
  float mat[order][order] = {{0}, {0}};
  uint8_t row, column;
  
  // Place matrix of type vector into an array
  for(row = 0; row < order; row++)
  {
    for(column = 0; column < order; column++)
    {
      mat[row][column] = matrix[row*order + column];
    }
  }

  //finding determinant
  uint8_t i = 0, j = 0, k = 0;
  k = 0;
  for(i = 0; i < order; i++)
  {
    determinant += mat[0][i]*(mat[(i+1)%order][(i+1)%order]*mat[(i+2)%order][(i+2)%order] - mat[(i+1)%order][(i+2)%order]*mat[(i+2)%order][(i+1)%order]);
  }

  if(!determinant)
  {
    Serial.println("Fail Inverse because determinant = 0");
    while(1){}
  }
  
  for(i = 0; i < order; i++)
  {
    for(j = 0; j < order; j++)
    {
			if(squared)
			{
				temp3x3[k] = ((mat[(j+1)%order][(i+1)%order] * mat[(j+2)%order][(i+2)%order]) - (mat[(j+1)%order][(i+2)%order] * mat[(j+2)%order][(i+1)%order]))/ determinant;
			}
			else
			{
				temp3[k] = ((mat[(j+1)%order][(i+1)%order] * mat[(j+2)%order][(i+2)%order]) - (mat[(j+1)%order][(i+2)%order] * mat[(j+2)%order][(i+1)%order]))/ determinant;
			}
      k++;
    }
  }
} */

void debugPrompt() {
	Serial.println(F("\n press 'g' and enter to continue\n"));
	char input;
	bool stopped = true;
	while(stopped)
	{
		if(Serial.available())
		{
			input = Serial.read();
			if(input == 'g')
			{
				stopped = false; 
			}
		}
	}
}

void print3x1(float matrix[]) {
	uint8_t i = 0;
	Serial.println();
	for(i = 0; i < 3; i++)
	{
		Serial.print(matrix[i], 6);
		Serial.print(F("\n"));
	}
	Serial.println();
}

void print3x3(float matrix[]) {
		uint8_t i = 0;
		Serial.println();
		for(i = 0; i < 9; i++)
		{
			Serial.print(matrix[i], 6);
			if(!((i+1)%3))
			{
				Serial.print(F("\n"));
			}
			else
			{
				Serial.print(F("\t"));
			}
		}
		Serial.println();
}

bool checkNaN3x1(float matrix[]) {
	uint8_t i = 0;
	for(i = 0; i < 3; i++)
	{
		if(isnan(matrix[i]))
		{
			return true;
		}
	}
	return false;
}

bool checkNaN3x3(float matrix[]) {
	uint8_t i = 0;
	for(i = 0; i < 9; i++)
	{
		if(isnan(matrix[i]))
		{
			return true;
		}
	}
	return false;
}

//	End of File
