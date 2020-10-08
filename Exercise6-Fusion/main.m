close all
clearvars
numMonteCarlo = 1000;
delT = 1; % Time step
F = kron([1 delT; 0 1],eye(2)); % State transition matrix
proNoise = 0.01; % Process noise
Q = proNoise*kron([(delT^3)/3 (delT^2)/2; (delT^2)/2 delT],eye(2)); %State transition noise matrix
sigmaRange = 5; % Measurement error standard deviation in range
sigmaTheta = 2*pi/180; % Measurement error standard deviation in angle
R = [sigmaRange^2 0;0 sigmaTheta^2]; % Measurement error covariance
radarState1 = [-600 800];
radarState2 = [600 800];
% For logging RMSE
ekfRmseTrackPos = zeros(1,60);
ekfRmseTrackVel = zeros(1,60);
ukfRmseTrackPos = zeros(1,60);
ukfRmseTrackVel = zeros(1,60);
rmseMeasPos1 = zeros(1,60);
rmseMeasPos2 = zeros(1,60);
% Loop through 500 Monte Carlo runs
for j = 1:numMonteCarlo
    % Generate a new target state
    targetState = zeros(4,60);
    targetState(:,1) = [0 0 0 10];
    for i = 2:60
        targetState(:,i) = F*targetState(:,i-1)+mvnrnd([0 0 0 0]',Q)';
    end
    measurements1 = zeros(2,60);
    measurements1Cart = zeros(2,60);
    measurements2 = zeros(2,60);
    measurements2Cart = zeros(2,60);
    ekfEstimate = zeros(4,60);
    ukfEstimate = zeros(4,60);
    for each time step
        for i = 1:60
            % Generate measurement from radar 1
            z1 = GenerateMeasurement(targetState(:,i),R,radarState1);
            measurements1(:,i) = z1;
            [measX1,measY1] = pol2cart(z1(2),z1(1));
            measurements1Cart(:,i) = [radarState1(1)+measX1;radarState1(2)+measY1];
            z2 = GenerateMeasurement(targetState(:,i),R,radarState2);
            measurements2(:,i) = z2;
            [measX2,measY2] = pol2cart(z2(2),z2(1));
            measurements2Cart(:,i) = [radarState2(1)+measX2;radarState2(2)+measY2];
            if i == 1
                ekfEstimate(:,i) = [0.5*(radarState1(1)+measX1)+0.5*(radarState2(1)+measX2); 0.5*(radarState1(2)+measY1)+0.5*(radarState2(2)+measY2);0;0];
                ukfEstimate(:,i) =[0.5*(radarState1(1)+measX1)+0.5*(radarState2(1)+measX2); 0.5*(radarState1(2)+measY1)+0.5*(radarState2(2)+measY2);0;0];
            elseif i == 2
            % Estimate initialisation
                ekfMean = [radarState1(1)+measX1;radarState1(2)+measY1;radarState1(1)+measX1- ekfEstimate(1,i-1);radarState1(2)+measY1-ekfEstimate(2,i-1)];
                ukfMean = [radarState1(1)+measX1;radarState1(2)+measY1;radarState1(1)+measX1-ukfEstimate(1,i-1);radarState1(2)+measY1-ukfEstimate(2,i-1)];
            % Covariance initialisation
                T = [cos(measurements1(2,i)) - measurements1(1,i)*sin(measurements1(2,i));sin(measurements1(2,i))measurements1(1,i)*cos(measurements1(2,i));];
                Rcart = T*R*T'; % Measurement error in Cartesian frame
                %ekfCovar = kron([1 1/delT; 1/delT 2/(delT*delT)],Rcart);
                %ukfCovar = kron([1 1/delT; 1/delT 2/(delT*delT)],Rcart);
                ekfCovar = [Rcart(1,1) Rcart(1,2) 0 0; Rcart(2,1)
                Rcart(2,2) 0 0; 0 0 2*Rcart(1,1) 0; 0 0 0 2*Rcart(2,2)];
                ukfCovar = [Rcart(1,1) Rcart(1,2) 0 0; Rcart(2,1)
                Rcart(2,2) 0 0; 0 0 2*Rcart(1,1) 0; 0 0 0 2*Rcart(2,2)];
                [ekfMean, ekfCovar] = EkfUpdate(ekfMean,ekfCovar,measurements2(:,i),R,radarState2);
                [ukfMean, ukfCovar] = UkfUpdate(ukfMean,ukfCovar,measurements2(:,i),R,radarState2);
                ekfEstimate(:,i) = ekfMean;
                ukfEstimate(:,i) = ukfMean;
            elseif i > 2
            % EKF Prediction
                ekfPriorMean = F*ekfMean;
                ekfPriorCovar = F*ekfCovar*F'+Q;
            % UKF Prediction
                ukfPriorMean = F*ukfMean;
                ukfPriorCovar = F*ukfCovar*F'+Q;
            % Extended Kalman update step
                [ekfMean, ekfCovar] = EkfUpdate(ekfPriorMean,ekfPriorCovar,measurements1(:,i),R,radarState1);
                [ekfMean, ekfCovar] = EkfUpdate(ekfMean,ekfCovar,measurements2(:,i),R,radarState2);
                ekfEstimate(:,i) = ekfMean;
                [ukfMean, ukfCovar] = UkfUpdate(ukfPriorMean,ukfPriorCovar,measurements1(:,i),R,radarState1);
                [ukfMean, ukfCovar] = UkfUpdate(ukfMean,ukfCovar,measurements2(:,i),R,radarState2);
                ukfEstimate(:,i) = ukfMean;
            end
        end
    % Log EKF RMSE
    ekfRmseTrackPos = ekfRmseTrackPos + (ekfEstimate(1,:)- targetState(1,:)).^2 + (ekfEstimate(2,:)-targetState(2,:)).^2;
    ekfRmseTrackVel = ekfRmseTrackVel + (ekfEstimate(3,:)- targetState(3,:)).^2 + (ekfEstimate(4,:)-targetState(4,:)).^2;
    % Log UKF RMSE
    ukfRmseTrackPos = ukfRmseTrackPos + (ukfEstimate(1,:)- targetState(1,:)).^2 + (ukfEstimate(2,:)-targetState(2,:)).^2;
    ukfRmseTrackVel = ukfRmseTrackVel + (ukfEstimate(3,:)- targetState(3,:)).^2 + (ukfEstimate(4,:)-targetState(4,:)).^2;
    % Log measurement RMSE
    rmseMeasPos1 = rmseMeasPos1 + (measurements1Cart(1,:)- targetState(1,:)).^2 + (measurements1Cart(2,:)-targetState(2,:)).^2;
    rmseMeasPos2 = rmseMeasPos2 + (measurements2Cart(1,:)- targetState(1,:)).^2 + (measurements2Cart(2,:)-targetState(2,:)).^2;
    if j==1
        figure
        plot(targetState(1,:),targetState(2,:))
        hold on
        plot(ekfEstimate(1,:),ekfEstimate(2,:),'r')
        plot(ukfEstimate(1,:),ukfEstimate(2,:),'g')
        plot(measurements1Cart(1,:),measurements1Cart(2,:),'x')
        plot(measurements2Cart(1,:),measurements2Cart(2,:),'kx')
        ylim([0 600])
        xlim([-200 200])
        xlabel('X (m)')
        ylabel('Y (m)')
        title('Target State, Estimate and Measurements')
    end
end
ekfRmseTrackPos = sqrt(ekfRmseTrackPos/numMonteCarlo);
ukfRmseTrackPos = sqrt(ukfRmseTrackPos/numMonteCarlo);
rmseMeasPos1 = sqrt(rmseMeasPos1/numMonteCarlo);
rmseMeasPos2 = sqrt(rmseMeasPos2/numMonteCarlo);
ekfRmseTrackVel = sqrt(ekfRmseTrackVel/numMonteCarlo);
ukfRmseTrackVel = sqrt(ukfRmseTrackVel/numMonteCarlo);
figure
plot(1:60,ekfRmseTrackPos(1:60),'r')
hold on
plot(1:60,ukfRmseTrackPos(1:60),'g')
plot(1:60,rmseMeasPos1(1:60))
plot(1:60,rmseMeasPos2(1:60),'k')
ylim([0 200])
xlabel('Time')
ylabel('RMSE (m)')
title('Root mean square error in position')
legend('EKF-Track','UKF-Track','Meas. Radar 1','Meas. Radar 2')
figure
plot(1:60,ekfRmseTrackPos(1:60),'r')
hold on
plot(1:60,ukfRmseTrackPos(1:60),'g')
plot(1:60,rmseMeasPos1(1:60))
plot(1:60,rmseMeasPos2(1:60),'k')
ylim([0 20])
xlabel('Time')
ylabel('RMSE (m)')
title('Root mean square error in position')
legend('EKF-Track','UKF-Track','Meas. Radar 1','Meas. Radar 2')
figure
plot(1:60,ekfRmseTrackVel(1:60),'r')
hold on
plot(1:60,ukfRmseTrackVel(1:60),'g')
ylim([0 50])
xlabel('Time')
ylabel('RMSE (m/s)')
title('Root mean square error in velocity')
legend('EKF-Track','UKF-Track')