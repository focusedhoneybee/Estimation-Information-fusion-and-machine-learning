function [ret]= MeasurementFunction(trueState,radarState)
    [h2,h1] = cart2pol(trueState(1)-radarState(1),trueState(2)- radarState(2));
    ret = [h1 h2]';
