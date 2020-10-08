function [xPosterior,PPosterior]= UkfUpdate(xPrior,PPrior,z,R,radarState)
    kappa = 0;
    na = length(xPrior);
    numberSigmaPoints = 1+2*na;
    sigmaPoints = zeros(na,numberSigmaPoints);
    transformedSigmaPoints = zeros(length(R),numberSigmaPoints);
    weights = zeros(1,numberSigmaPoints);
    sigmaPoints(:,1) = xPrior;
    weights(1) = kappa/(na+kappa);
    for i = 1:na
        root = sqrtm((na+kappa)*PPrior);
        sigmaPoints(:,i+1) = xPrior + root(i,:)';
        weights(i+1) = 1/(2*(na+kappa));
    end
    for i = 1:na
        root = sqrtm((na+kappa)*PPrior);
        sigmaPoints(:,na+1+i) = xPrior - root(i,:)';
        weights(na+1+i) = 1/(2*(na+kappa));
    end
    for i = 1:numberSigmaPoints
        transformedSigmaPoints(:,i) = MeasurementFunction(sigmaPoints(:,i),radarState);
    end
    zHat = zeros(length(R),1);
    for i = 1:numberSigmaPoints
        zHat = zHat + weights(i)*transformedSigmaPoints(:,i);
    end
    Pxz = zeros(na,length(R));
    for i = 1:numberSigmaPoints
        Pxz = Pxz + weights(i)*(sigmaPoints(:,i)-
        xPrior)*(transformedSigmaPoints(:,i)-zHat)';
    end
    Pzz = zeros(length(R),length(R));
    for i = 1:numberSigmaPoints
        Pzz = Pzz + weights(i)*(transformedSigmaPoints(:,i)- zHat)*(transformedSigmaPoints(:,i)-zHat)';
    end
    S = R + Pzz;
    K = Pxz/S;
    xPosterior = xPrior + K*(z-zHat);
    PPosterior = PPrior - K*S*K';