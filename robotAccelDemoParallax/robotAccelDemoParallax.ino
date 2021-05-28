/*
 * Author:        Jacob Morgan
 * Date Created:  09/25/2020
 * Company:       UNCC
 * Extra:         Using modified library of Adafruit MMA8451 library
 * 
 * Description:   Using Adafruit MMA8451 library, acceleromters are read
 *                to determine wheel position. An extended kalman filter is run to
 *                significantly filter the data for usable motion estimates.
*/

#include <TimedAction.h>
#include <Adafruit_MMA8451.h>
#include <Adafruit_Sensor.h>
#include <Kalman.h>
#include <Servo.h>
#include <WString.h>

#define LEFTSERVOPIN  13
#define RIGHTSERVOPIN 12

#define STD_DEV_SENSOR      0.87
#define STD_DEV_MODEL       0.0294935

#define SENSORWINDOWSIZE    5
#define TIMEWINDOW          100

#define SENSORTIME_ms       0002        // Time interval to check sensors
#define PRINTTIME_ms        0004

#define MICROS_TO_SEC       0.000001

#define WHEELDIAMETER       0.066       // diameter of wheel with tire in meters

#define WHEELCIRCUMFERENCE PI*WHEELDIAMETER

//#define DEBUG

// function prototype
void stopWheels();
void moveForward();
//void moveBackward();
//void turnLeft();
//void turnRight();
void averageInitialValue();
void getSensorValues();
void displayValues();
void averageDeltaTime();

// Kalman and Acceleromter objects
Kalman            leftSensor  = Kalman();
Adafruit_MMA8451  leftAccel   = Adafruit_MMA8451();
sensors_event_t   eventL;

unsigned long startMicros     = 0;
unsigned long currentTime     = 0;
unsigned long leftTime        = 0;                    // Number of ticks since last left sensor read
unsigned int  intLeftTime     = 0;

float leftDeltaT              = 0;                    // Ticks converted to time value
float leftAccelValueX         = 0;
float leftAccelValueY         = 0;
float leftWindowX[SENSORWINDOWSIZE]   = {0};
float leftWindowY[SENSORWINDOWSIZE]   = {0};
float timeWindow[TIMEWINDOW] = {0};

uint8_t timeIndex = 0;

char input;
bool userStop;
bool motorRunning;

Servo leftServo;
Servo rightServo;

TimedAction sensorAction  = TimedAction(SENSORTIME_ms, getSensorValues);
TimedAction printAction   = TimedAction(PRINTTIME_ms, displayValues);

//----------------------------------------------------------------------
//------------Set up of Arduino-----------------------------------------
//----------------------------------------------------------------------
void setup(void) 
{
  leftServo.attach(LEFTSERVOPIN);
  rightServo.attach(RIGHTSERVOPIN);
  stopWheels();
  Serial.begin(115200);

  if (!leftAccel.begin(MMA8451_LEFTWHEEL))
  {
    Serial.println(F("Couldnt start left accelerometer"));
    while (1){}
  }

  // Init the accelerometer
  leftAccel.setRange(MMA8451_RANGE_2_G);
  leftAccel.setDataRate(MMA8451_DATARATE_800_HZ);

  // Set values for Kalman Objects
  // Set more trust on filter by providing lower standard deviation value
  leftSensor.initialize();
  leftSensor.setStdDevModel(STD_DEV_MODEL);
  leftSensor.setStdDevSensor(STD_DEV_SENSOR);
  leftSensor.setWheelCircumference_m(WHEELCIRCUMFERENCE);
  leftSensor.setSinusoidBounds(12.5);
} // end of setup

//----------------------------------------------------------------------
//------------Main Forever Loop-----------------------------------------
//----------------------------------------------------------------------
void loop()
{
  userStop = true;
  Serial.println(F("'w' = forward, 's' = Backward,")); // User prompt
  Serial.println(F("'a' = left,    'd' = right")); // User prompt
  while(userStop)
  {
    if(Serial.available())
    {
      input = Serial.read();
      if(input == 'w' || input == 's' || input == 'a' || input == 'd')
      {
        userStop = false; 
      }
    }
  }
  Serial.println(F("enter 's' + send to stop"));
  Serial.println(F("\n\n"));
  Serial.println(F("us\tLx\tkx"));     // Print column titles (time, accel value, kalman value)
  Serial.flush();
  delay(1000);

//--------------------------------------------------------
//------------Movement Portion of loop--------------------
//--------------------------------------------------------
  averageInitialValue();      // Average inital x and y accel values
  startMicros = micros();     // Get initial start time for motors
  leftTime  = startMicros;
  Serial.print(F("Inital Average X =\t")); Serial.println(leftAccelValueX, 6);

  switch(input)
  {
    case 'w':
      leftSensor.setWheelDirection(Kalman::direction::COUNTERCLOCK);
      moveForward();
      break;/*
    case 's':
      leftSensor.setWheelDirection(Kalman::direction::CLOCKWISE);
      moveBackward();
      break;
    case 'a':
      leftSensor.setWheelDirection(Kalman::direction::CLOCKWISE);
      turnLeft();
      break;
    case 'd':
      leftSensor.setWheelDirection(Kalman::direction::COUNTERCLOCK);
      turnRight();
      break;*/
    default:
      break;
  }
  input = 0;

  leftSensor.initKalmanPrediction(leftAccelValueX, leftAccelValueY);

  sensorAction.enable();
  printAction.enable();
  motorRunning = true;
  while(motorRunning)
  {
    if(Serial.available())      // Check to stop motors
    {
      if(Serial.read() == 's')
      {
        motorRunning = false; 
      }
    }

    if(leftSensor.getAbortStatus())
    {
      stopWheels();
      motorRunning = false;
    }

    sensorAction.check();       // Check to read sensor values and run kalman filters
    printAction.check();        // Check to print out values
  }
  stopWheels();
  Serial.println();             // Add line space between next user prompt
//  averageDeltaTime();
  Serial.flush();
  while(1){}
} // end forever loop

//----------------------------------------------------------------------
//------------Display Sensor and Filtered Data--------------------------
//----------------------------------------------------------------------
void displayValues()    // Called when printAction times out
{
  Serial.print(intLeftTime);              Serial.print(F("\t"));
  Serial.print(eventL.acceleration.x, 6); Serial.print(F("\t"));
  Serial.println(leftSensor.getPredictionValueX(), 6);
}

//----------------------------------------------------------------------
//------------Read Sensor Values of Accelerometers----------------------
//----------------------------------------------------------------------
void getSensorValues()  // Called when sensorAction times out
{
  leftAccel.read();
  currentTime = micros();
  intLeftTime = currentTime - leftTime;               // Current time - last left time
  leftDeltaT  = ((float) intLeftTime)*MICROS_TO_SEC;  // Get time difference as float in seconds
  leftTime    = currentTime;                          // Set last left time to current time

  leftAccel.getEvent(&eventL);
  leftAccelValueX = eventL.acceleration.x;

  timeWindow[timeIndex] = leftDeltaT;
  timeIndex++;
  timeIndex%=TIMEWINDOW;
      
  // set sensor values in Kalman filter and run
  leftSensor.setSensorReadTimeDelta(leftDeltaT);
  leftSensor.setSensorReadings(leftAccelValueX);
  leftSensor.runFilter();
} // end get sensor values

//----------------------------------------------------------------------
//------------Stop Both Wheel-------------------------------------------
//----------------------------------------------------------------------
void stopWheels()
{
  leftServo.writeMicroseconds(1520);  // Stop Wheels
  rightServo.writeMicroseconds(1520);

  sensorAction.disable();             // Turn off timers for running print and sensor reads
  printAction.disable();

  leftSensor.resetKalmanFilter();     // Reset Kalman values
}

//----------------------------------------------------------------------
//------------Go Straight-----------------------------------------------
//----------------------------------------------------------------------
void moveForward()
{
  leftSensor.setWheelPeriod_s(1.565);
  leftSensor.setMaxMinSensorValues(10.299857, 10.433932);
  leftServo.writeMicroseconds(1650);
  rightServo.writeMicroseconds(1390);
}

/*
//----------------------------------------------------------------------
//------------Go Backward-----------------------------------------------
//----------------------------------------------------------------------
void moveBackward()
{
  leftServo.writeMicroseconds(1440);
  leftSensor.setWheelPeriod_s(1.7275);
}

//----------------------------------------------------------------------
//------------Turn Left-------------------------------------------------
//----------------------------------------------------------------------
void turnLeft()
{
  leftServo.writeMicroseconds(1310);
  leftSensor.setWheelPeriod_s(1.5764);
}

//----------------------------------------------------------------------
//------------Turn Right------------------------------------------------
//----------------------------------------------------------------------
void turnRight()
{
  leftServo.writeMicroseconds(1700);
  leftSensor.setWheelPeriod_s(1.5764);
}
*/
//----------------------------------------------------------------------
//------------Average initial value-------------------------------------
//----------------------------------------------------------------------
void averageInitialValue()
{
  leftAccelValueX = 0;
  leftAccelValueY = 0;
  
  for(int i = 0; i < SENSORWINDOWSIZE; i++)
  {
    leftAccel.read();
    leftAccel.getEvent(&eventL);
    leftAccelValueX += eventL.acceleration.x;
    leftAccelValueY += eventL.acceleration.y;
    delay(2);
  }

  leftAccelValueX /= SENSORWINDOWSIZE;
  leftAccelValueY /= SENSORWINDOWSIZE;
}

//----------------------------------------------------------------------
//------------Average deltaT value--------------------------------------
//----------------------------------------------------------------------
void averageDeltaTime()
{
  double sum = 0;
  for(int i = 0; i < TIMEWINDOW; i++)
  {
    sum += (double) timeWindow[i];
  }
  sum /= TIMEWINDOW;
  Serial.print(F("\nAverage time window =\t")); Serial.println(sum, 6);
}
