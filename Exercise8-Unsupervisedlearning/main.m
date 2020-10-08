close all
clear all
load data.mat
% To store the current assignments
assignments = zeros(1,length(data));
% 'Centroids' stores the current centroid
centroidIndex = randi(length(data),2,1);
centroids = data(:,centroidIndex);
while 1
    %Update assignments
    for j = 1:length(data)
        squaredMagnitude = sum((centroids -repmat(data(:,j),1,2)).^2);
        [y, index] = min(squaredMagnitude);
        assignments(j) = index;
    end
    lastCentroids = centroids;
    %Update centroid
    for j = 1:length(data)
        for k = 1:2
        samples = data(:,assignments==k);
        centroids(:,k) = mean(samples,2);
        end
    end
    % Exit the loop on convergence
    if lastCentroids == centroids
        break;
    end
end
figure
plot(data(1,assignments==1),data(2,assignments==1),'bo')
hold on
plot(data(1,assignments==2),data(2,assignments==2),'rx')
xlabel('x1')
ylabel('x2')
title('Result of K-Means Clustering')


% For storing the posterior probabilities
posteriorProbs = zeros(2,length(data));
% For storing the model parameters
phi = [0.5;0.5];
meansIndex = randi(length(data),2,1);
means = data(:,meansIndex);
covariances = zeros(2,2,2);
covariances(:,:,1) = [1 0; 0 1];
covariances(:,:,2) = [1 0; 0 1];
while 1
    % Update the posterior probabilities
    for j = 1:length(data)
        for m = 1:2
            posteriorProbs(m,j) =
            mvnpdf(data(:,j),means(:,m),covariances(:,:,m))*phi(m);
        end
        posteriorProbs(:,j) = posteriorProbs(:,j)/
        sum(posteriorProbs(:,j));
    end
    % Update the model parameters
    lastPhi = phi;
    phi = sum(posteriorProbs,2) / length(data);
    means(:,1) = sum(posteriorProbs(1,:).*data(:,:),2)/sum(posteriorProbs(1,:));
    means(:,2) = sum(posteriorProbs(2,:).*data(:,:),2)/sum(posteriorProbs(2,:));
    covariances = zeros(2,2,2);
    for j=1:length(data)
        covariances(:,:,1) = covariances(:,:,1) +
        posteriorProbs(1,j)*(data(:,j) - means(:,1)) * (data(:,j) -means(:,1))';
        covariances(:,:,2) = covariances(:,:,2) +posteriorProbs(2,j)*(data(:,j) - means(:,2)) * (data(:,j) -means(:,2))';
    end
    covariances(:,:,1) = covariances(:,:,1) /
    sum(posteriorProbs(1,:));
    covariances(:,:,2) = covariances(:,:,2) /
    sum(posteriorProbs(2,:));
    if lastPhi == phi
        break;
    end
end
