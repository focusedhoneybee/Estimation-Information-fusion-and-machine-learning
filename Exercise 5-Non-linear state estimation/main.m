close all
clearvars

% Set the number of Monte Carlo runs
numMonteCarlo = 1000;

delT = 1; % Time step
% Task 1 - Complete the state transition matrix
F = []; % State transition matrix
proNoise = 0.01; % Process noise
% Task 1 - Complete the state transition noise matrix
Q = []; % State transition noise matrix

sigmaRange = 5; % Measurement error standard deviation in range
sigmaTheta = 4*pi/180; % Measurement error standard deviation in angle
% Task 2 - Complete the measurement error covariance
R = []; % Measurement error covariance

radarState = [-600 800];

% Preallocate arrays for holding RMSE data
ekfRmseTrackPos = zeros(1,60);
ekfRmseTrackVel = zeros(1,60);
ukfRmseTrackPos = zeros(1,60);
ukfRmseTrackVel = zeros(1,60);
rmseMeasPos = zeros(1,60);

% Loop through the Monte Carlo runs
for j = 1:numMonteCarlo
  
    % Generate a random target trajectory according to the model
    targetState = zeros(4,60);
    targetState(:,1) = [0 0 0 10]; % Starting state
    for i = 2:60
        % Perform a random walk according to motion model
        targetState(:,i) = F*targetState(:,i-1)+mvnrnd([0 0 0 0]',Q)'; 
    end
    
    % Arrays for logging data from a single run
    measurements = zeros(2,60); % Log of measurements in polar space
    measurementsCart = zeros(2,60); % Log of the measurements in Cartesian space
    ekfEstimate = zeros(4,60); % Log of EKF estimate
    ukfEstimate = zeros(4,60); % Log of UKF estimate
    
    % Loop for each time step
    for i = 1:60
        % Generate a random measurement
        z = GenerateMeasurement(targetState(:,i),R,radarState);   
        % Store the measurement
        measurements(:,i) = z; 
        [measX,measY] = pol2cart(z(2),z(1)); % Convert measurement into Cartesian space       
        measurementsCart(:,i) = [radarState(1)+measX;radarState(2)+measY]; % Store measurement
        if i == 1        
            % Store the first time step to be used later in two point
            % initialisation
            ekfEstimate(:,i) = [radarState(1)+measX; radarState(2)+measY;0;0];
            ukfEstimate(:,i) = [radarState(1)+measX; radarState(2)+measY;0;0];
        elseif i == 2
            % Perform two point initialisation
            ekfMean = [radarState(1)+measX;radarState(2)+measY;radarState(1)+measX-ekfEstimate(1,i-1);radarState(2)+measY-ekfEstimate(2,i-1)];
            ukfMean = [radarState(1)+measX;radarState(2)+measY;radarState(1)+measX-ukfEstimate(1,i-1);radarState(2)+measY-ukfEstimate(2,i-1)];
            ekfEstimate(:,i) = ekfMean;            
            ukfEstimate(:,i) = ukfMean;            
            % Covariance initialisation
            T = [
                cos(measurements(2,i)) -measurements(1,i)*sin(measurements(2,i));
                sin(measurements(2,i)) measurements(1,i)*cos(measurements(2,i));
                ];
            Rcart = T*R*T'; % Measurement error in Cartesian frame
            ekfCovar = [Rcart(1,1) Rcart(1,2) 0 0; Rcart(2,1) Rcart(2,2) 0 0; 0 0 2*Rcart(1,1) 0; 0 0 0 2*Rcart(2,2)];
            ukfCovar = [Rcart(1,1) Rcart(1,2) 0 0; Rcart(2,1) Rcart(2,2) 0 0; 0 0 2*Rcart(1,1) 0; 0 0 0 2*Rcart(2,2)];
        elseif i > 2
            % EKF Prediction
            [ekfPriorMean,ekfPriorCovar] = kalmanPrediction(ekfMean,ekfCovar,F,Q);
            % UKF Prediction
            [ukfPriorMean,ukfPriorCovar] = kalmanPrediction(ukfMean,ukfCovar,F,Q);
            % Extended Kalman update step
            [ekfMean, ekfCovar] = EkfUpdate(ekfPriorMean,ekfPriorCovar,z,R,radarState);
            ekfEstimate(:,i) = ekfMean;
            % Extended Kalman update step
            [ukfMean, ukfCovar] = UkfUpdate(ukfPriorMean,ukfPriorCovar,z,R,radarState);
            ukfEstimate(:,i) = ukfMean;
        end
    end
    % Log EKF RMSE
    ekfRmseTrackPos = ekfRmseTrackPos + (ekfEstimate(1,:)-targetState(1,:)).^2 + (ekfEstimate(2,:)-targetState(2,:)).^2;
    ekfRmseTrackVel = ekfRmseTrackVel + (ekfEstimate(3,:)-targetState(3,:)).^2 + (ekfEstimate(4,:)-targetState(4,:)).^2;
    % Log UKF RMSE
    ukfRmseTrackPos = ukfRmseTrackPos + (ukfEstimate(1,:)-targetState(1,:)).^2 + (ukfEstimate(2,:)-targetState(2,:)).^2;
    ukfRmseTrackVel = ukfRmseTrackVel + (ukfEstimate(3,:)-targetState(3,:)).^2 + (ukfEstimate(4,:)-targetState(4,:)).^2;
    % Log measurement RMSE
    rmseMeasPos = rmseMeasPos + (measurementsCart(1,:)-targetState(1,:)).^2 + (measurementsCart(2,:)-targetState(2,:)).^2;
    % Make a plot of a single run for only the single run
    if j==1
        figure
        plot(targetState(1,:),targetState(2,:))
        hold on
        plot(ekfEstimate(1,:),ekfEstimate(2,:),'r')
        plot(ukfEstimate(1,:),ukfEstimate(2,:),'g')
        plot(measurementsCart(1,:),measurementsCart(2,:),'x')
        ylim([0 600])
        xlim([-200 200])
        xlabel('X (m)')
        ylabel('Y (m)')
        title('Target State, Estimate and Measurements')
    end
end

% Calculate the RMSE
ekfRmseTrackPos = sqrt(ekfRmseTrackPos/numMonteCarlo);
ukfRmseTrackPos = sqrt(ukfRmseTrackPos/numMonteCarlo);
rmseMeasPos = sqrt(rmseMeasPos/numMonteCarlo);
ekfRmseTrackVel = sqrt(ekfRmseTrackVel/numMonteCarlo);
ukfRmseTrackVel = sqrt(ukfRmseTrackVel/numMonteCarlo);

% Plot the RMSE in position
figure
plot(1:60,ekfRmseTrackPos(1:60),'r')
hold on
plot(1:60,ukfRmseTrackPos(1:60),'g')
plot(1:60,rmseMeasPos(1:60))
ylim([0 200])
xlabel('Time')
ylabel('RMSE (m)')
title('Root mean square error in position')
legend('EKF-Track','UKF-Track','Measurements')

% Plot the RMSE in velocity
figure
plot(1:60,ekfRmseTrackVel(1:60),'r')
hold on
plot(1:60,ukfRmseTrackVel(1:60),'g')
ylim([0 50])
xlabel('Time')
ylabel('RMSE (m/s)')
title('Root mean square error in velocity')
legend('EKF-Track','UKF-Track')