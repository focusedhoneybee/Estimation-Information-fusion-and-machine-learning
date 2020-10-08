T = 1; % Time step is one second
F = [1 0 T 0; 0 1 0 T; 0 0 1 0; 0 0 0 1]; % State transition matrix
proNoise = 0.01; % Process noise intensity q
Q = proNoise*[(T^3)/3 0 (T^2)/2 0 ; 0 (T^3)/3 0 (T^2)/2; (T^2)/2 0 T
 0; 0 (T^2)/2 0 T]; % Noise matrix
sigmaX = 5; % Measurement error standard deviation in x
sigmaY = 5; % Measurement error standard deviation in y
R = [sigmaX^2 0;0 sigmaY^2]; % Measurement error covariance matrix
H = [1 0 0 0; 0 1 0 0]; % Observation matrix
% Arrays for storing a deterministic and four non-deterministic
% trajectories
targetStateDeter = zeros(4,60);
targetStateDeter(:,1) = [0 0 0 10];
targetState = zeros(4,60);
targetState(:,1) = [0 0 0 10];
targetState2 = zeros(4,60);
targetState2(:,1) = [0 0 0 10];
targetState3 = zeros(4,60);
targetState3(:,1) = [0 0 0 10];
targetState4 = zeros(4,60);
targetState4(:,1) = [0 0 0 10];
% Generate trajectories
for i = 2:60
 % Deterministic (without noise)
 targetStateDeter(:,i) = F*targetStateDeter(:,i-1);
 % Non-deterministic (with noise)
 targetState(:,i) = F*targetState(:,i-1)+mvnrnd([0 0 0 0]',Q)';
 targetState2(:,i) = F*targetState2(:,i-1)+mvnrnd([0 0 0 0]',Q)';
 targetState3(:,i) = F*targetState3(:,i-1)+mvnrnd([0 0 0 0]',Q)';
 targetState4(:,i) = F*targetState4(:,i-1)+mvnrnd([0 0 0 0]',Q)';
end
% Plot for exercise 1
figure
plot(targetStateDeter(1,:),targetStateDeter(2,:))
xlabel('X (m)')
ylabel('Y (m)')
title('Deterministic trajectory')
% Plot for exercise 2
figure
plot(targetStateDeter(1,:),targetStateDeter(2,:))
hold on
plot(targetState(1,:),targetState(2,:),'k')
plot(targetState2(1,:),targetState2(2,:),'k')
plot(targetState3(1,:),targetState3(2,:),'k')
plot(targetState4(1,:),targetState(2,:),'k')
ylim([0 600])
xlim([-60 60])
xlabel('X (m)')
ylabel('Y (m)')
title('Non-Deterministic trajectory')
% Arrays for storing measurements, estimate from track and Kalman
 filter
% gain
measurements = zeros(2,60);
% Indexing for 60 times steps
for i = 1:60
 % Generate a measurement based on
 z = H*targetState(:,i)+mvnrnd([0 0]',R)';
 %Store all the measurements
 measurements(:,i) = z;
end
% Plot for exercise 3
figure
plot(targetState(1,:),targetState(2,:))
hold on
plot(measurements(1,:),measurements(2,:),'x')
ylim([0 600])
xlim([-60 60])
xlabel('X (m)')
ylabel('Y (m)')
title('Target State and Measurements')