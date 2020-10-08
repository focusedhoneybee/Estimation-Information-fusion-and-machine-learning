function [xPosterior,PPosterior]= EkfUpdate(xPrior,PPrior,z,R,radarState)
    H = MeasurementMatrix(xPrior,radarState);
    zHat = MeasurementFunction(xPrior,radarState);
    S = H*PPrior*H'+R;
    v = z - zHat;
    W = PPrior*H'*pinv(S);
    xPosterior = xPrior+W*v;
    PPosterior = PPrior-W*S*W';
end