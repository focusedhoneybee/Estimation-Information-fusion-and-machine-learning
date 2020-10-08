close all
clearvars

%trueRange = 5.25;
%trueMeanSnr = measurementFunction(trueRange);
%data = -trueMeanSnr*log(rand(20,1));

load data.mat
parameterRange = 0:0.01:10;

% Plot the likelihood as a function of range
figure
plot(parameterRange,likelihood(data,parameterRange));
xlabel('Range (m)')
ylabel('Likelihood')

% Find maximum likelihood estimate
funcHandle = @(r)-likelihood(data,r);
mleEstimate = fminbnd(funcHandle,0,10);
disp(['ML Estimate: ', num2str(mleEstimate)])

% Plot the posterior as a function of range
priorMean = 5;
priorStd = 1;
funcHandle = @(r)likelihood(data,r).*normpdf(r,priorMean,priorStd);
posteriorDenominator = integral(funcHandle,0,Inf);
figure
posterior = likelihood(data,parameterRange).*normpdf(parameterRange,priorMean,priorStd)/
posteriorDenominator;
plot(parameterRange,posterior);
xlabel('Range (m)')
ylabel('Probability Density')

% Find maximum as posteriori estimate
funcHandle = @(r)-likelihood(data,r).*normpdf(r,priorMean,priorStd);
mapEstimate = fminbnd(funcHandle,0,10);
disp(['MAP Estimate: ', num2str(mapEstimate)])

% Find least squares estimate
funcHandle = @(r)sum((data - measurementFunction(r)).^2);
lsEstimate = fminbnd(funcHandle,0,10);
disp(['LS Estimate: ', num2str(lsEstimate)])

% Find the MMSE estimate
funcHandle = @(r)r.*likelihood(data,r).*normpdf(r,priorMean,priorStd)/
posteriorDenominator;
mmseEstimate = integral(funcHandle,0,Inf);
disp(['MMSE Estimate: ', num2str(mmseEstimate)])