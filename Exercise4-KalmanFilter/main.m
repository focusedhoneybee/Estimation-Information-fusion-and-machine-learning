close all
clearvars

T = 1; % Time step is one second
F = []; % TASK 1 - Complete the state transition matrix
proNoise = 0.01; % Process noise intensity q
Q = proNoise*[]; % TASK 1 - Complete the transition noise covariance matrix

sigmaX = 5; % Measurement error standard deviation in x
sigmaY = 5; % Measurement error standard deviation in y
R = []; % TASK 2 - Complete the measurement error covariance matrix
H = []; % TASK 2 - Complete the measurement matrix

load('data.mat')

estimate = zeros(4,60);
gain = zeros(1,60);

% Indexing for 60 times steps
for i = 1:60 
   %Store all the measurements
   z = measurements(:,i); 
   if i == 2
      % In the second time step perform initialisation
      mean = [z(1) z(2) z(1)-measurements(1,i-1) z(2)-measurements(2,i-1)]';
      covar = [R(1,1) 0 R(1,1) 0; 0 R(2,2) 0 R(2,2); R(1,1) 0 2*R(1,1) 0; 0 R(2,2) 0 2*R(2,2)];
      estimate(:,i) = mean;
   elseif i > 2
       
       % Perform the Kalman filter prediction
       [priorMean, priorCovar] = kalmanPrediction(mean,covar,F,Q);
       
       if i == 4
           % TASK 4 - Plot the prior pdf using surf and mvnpdf        
       end
       
       % Perform the Kalman filter update and log the Kalman gain
       % additionally
       [mean,covar,gain(:,i)] = kalmanUpdate(priorMean,priorCovar,z,H,R);
       
       if i == 4
           % TASK 4 - Plot the posterior pdf using surf and mvnpdf  
       end

       % Log the estimate
       estimate(:,i) = mean;
   end
end

% TASK 5 - Plot the true state, the measurements and state estimate
figure
plot(targetState(1,:),targetState(2,:))
hold on
plot(estimate(1,:),estimate(2,:),'r')
plot(measurements(1,:),measurements(2,:),'x')
ylim([0 600])
xlim([-100 100])
xlabel('X (m)')
ylabel('Y (m)')
title('Target State, Estimate and Measurements')

% TASK 6 - Plot the Kalman filter gain
figure
plot(3:60,gain(3:60))
xlabel('Time (sec)')
ylabel('Kalman Gain')
title('Kalman filter gain')

