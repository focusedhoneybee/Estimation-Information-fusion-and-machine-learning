function [z]= GenerateMeasurement(trueState,R,radarState)

noise = mvnrnd([0 0]',R)';
z = MeasurementFunction(trueState,radarState) + noise;