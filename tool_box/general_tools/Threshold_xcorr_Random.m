% This function assign random positions to time2 vector and construct a ccg
% using the times1 vector. The output is a distribution of maximal values
% coming from this procedure.
%
% output = Threshold_xcorr_Random(times1, times2, binSize, duration, smooth, mode, iterations)
%
% --- Inputs ---
% singla1: vector, data1
% signal2: vector, data2, they will be shuffle
% lags: int, number of total lags to calculate the xcorr
% mode: Str, 'coeff' or 'none', see xcorr 'mode' from Matlab internal f(x)
% iterations: int, times the ccg with shuffle data will be performed
%
% see xcorr from Matlab signal processing package
%
% --- OUTPUT ---
% output: vector containing the maximal value at each iteration.
% Morici Juan Facundo, 07/2023

function output = Threshold_xcorr_Random(signal1, signal2, lags, mode, iterations)

output = zeros(iterations,1);
for S = 1:iterations

    x = signal1;
    y = signal2;
    y = signal2(randperm(max(size(signal2))));
    
    if strcmp(mode,'coeff')
        c = xcorr(x,y,lags,'coeff');
        output(S) = max(c); clear c
    elseif strcmp(mode,'none')
        c = xcorr(x,y,lags);
        output(S) = max(c); clear c        
    end
end
end