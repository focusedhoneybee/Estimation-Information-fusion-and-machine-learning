function [priorMean,priorCovariance] = kalmanPrediction(posteriorMean,posteriorCovar,F,Q)
    priorMean = F*posteriorMean;
    priorCovariance = F*posteriorCovar*F'+Q;
end