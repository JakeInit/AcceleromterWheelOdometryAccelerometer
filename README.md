# AcceleromterWheelOdometryAccelerometer
The Arduino code for setting and running an Extended Kalman Filter when placing an accelerometer on a wheel. This code focuses on single axis monitoring and estimation.

How to run:
In Setup:
1. initialize()
2. setWheelCircumference_m(float wheelCircumference_m_);
3. setStdDevModel(float stdDev);
4. setStdDevSensor(float stdDev);
5. setSinusoidBounds(float bound);

When turning on motor to turn wheel
1. setWheelPeriod_s(float period_s);
2. setWheelDirection(direction wheelDirection_);
3. setMaxMinSensorValues(float min_, float max_);
4. initKalmanPrediction(float predictionX, float predictionY);

Then can turn on motors after doing above steps and entering the loop
Then in loop:
1. setSensorReadTimeDelta(float timeDelta_s);
2. setSensorReadings(float sensorReading_mpss);
3. runFilter();
3. get Kalman Results using getPredictionValueX() and getPredictionValueY()
5. get distance traveled using getDistanceTaveled_m()
