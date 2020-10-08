function [posteriorMean,posteriorCovariance,gainX] = kalmanUpdate(priorMean,priorCovariance,z,H,R)
    % Expected measurement mean
    z_hat = H*priorMean;
    % Expected measurement covariance
    S = R + H*priorCovariance*H';
    % Measurement residual
    v = z - z_hat;
    % Kalman filter weight
    W = priorCovariance*H'*inv(S);
    % Output the gain in x
    gainX = W(1,1);
    % New estimate mean
    posteriorMean = priorMean+W*v;
    % New estimate covariance
    posteriorCovariance = priorCovariance-W*S*W';
end