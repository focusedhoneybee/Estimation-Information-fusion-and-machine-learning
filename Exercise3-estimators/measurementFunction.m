function meanSnr = measurementFunction(range)
    nominalRange = 10;
    meanSnr = (nominalRange./range).^2;
end