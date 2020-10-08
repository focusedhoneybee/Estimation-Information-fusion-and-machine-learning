function [outputArg] = likelihood(snrSamples,range)
    outputArg = zeros(1,length(range));
    for i = 1:length(range)
        meanSnr = measurementFunction(range(i));
        outputArg(i) = prod((1/meanSnr)*exp(-snrSamples/meanSnr));
    end
end