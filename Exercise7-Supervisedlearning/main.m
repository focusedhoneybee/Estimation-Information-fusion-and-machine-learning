close all
clear all
% Training data
x = [1 8 2;
 1 16 16;
 1 3.2 12;
 1 20 6;
 1 14 10;
 1 10 8];
% Labels
y =[0.0 1.0 0.0 1.0 1.0 0.0];
alpha = 0.1; % Training rate
% Plot the data
figure
hold on
for i = 1:size(x,1)
    if(y(i) == 0.0)
        plot(x(i,2),x(i,3),'bo')
    else
        plot(x(i,2),x(i,3),'rx')
    end
end
xlim([0,20])
ylim([0,20])
xlabel('Car Age (Years)')
ylabel('Kilometers Driven (10000km)')
theta = [1 1 1]'; % Initial weight vector
while 1
    oldTheta = theta;
    for j = 1:size(x,1)
        if (theta'*x(j,:)' >= 0)
            modelOut = 1.0;
            else
            modelOut = 0.0;
        end
        if(modelOut ~= y(j))
            theta = theta +alpha*(y(j) - modelOut)*x(j,:)';
        end
    end
    if(oldTheta == theta)
        break
    end
end
% Plot the decision boundary for PLA
x1Vals = 0:0.1:20;
x2Vals = (-theta(1) -theta(2)*x1Vals)/theta(3);
plot(x1Vals,x2Vals)
