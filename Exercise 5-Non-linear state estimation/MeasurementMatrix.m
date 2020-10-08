function [measMatrix]= MeasurementMatrix(state,radarState)
    h11 = (state(1)-radarState(1))/((state(1)-radarState(1))^2 + (state(3)-radarState(2))^2)^0.5;
    h12 = (state(3)-radarState(2))/((state(1)-radarState(1))^2 + (state(3)-radarState(2))^2)^0.5;
    h21 = -(state(3)-radarState(2))/((state(1)-radarState(1))^2 + (state(3)-radarState(2))^2);
    h22 = (state(1)-radarState(1))/((state(1)-radarState(1))^2 + (state(3)-radarState(2))^2);
    measMatrix = [h11 h12 0 0; h21 h22 0 0];
