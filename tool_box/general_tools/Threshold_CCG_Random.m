% This function asign random positions to time2 vector and construct a ccg
% using the times1 vector. The output is a distribution of maximal values
% coming from this procedure.
%
% output = Threshold_CCG_Random(times1, times2, binSize, duration, smooth, mode, iterations)
%
% --- Inputs ---
% times1: vector, event1 times
% times2: vector, event2 times, they will be shuffle
% binSize: int/float, bin duration for CCG in sec
% duration: float/int, total duration of the CCG
% smooth: int for smoothing the CCG
% mode: Str, 'ccg' or 'ccv', see ccg funtion from FMAToolbox
%       if is 'ccg' the values are in probability
%       if is 'ccv' is a.u
% iterations: int, times the ccg with shuffle data will be performed
% Xedges: int for spatial bins for the firing curve
%
% --- OUTPUT ---
% output: vector containing the maximal value at each iteration.
% Morici Juan Facundo, 07/2023

function output = Threshold_CCG_Random(times1, times2, binSize, duration, smooth, mode, iterations)

output = zeros(iterations,1);
for S = 1:iterations
    a = min(times2) + rand*10;% 1° tiempo
    isi = diff(times2);
    isi_per = isi(randperm(length(isi)));%shuffle isi
    times2 = [a; isi_per];
    times2 =cumsum(times2);    
    
    x = times1;
    y = times2;
    [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
    [ccg,tttt] = CCG(s,ids,'binSize',binSize,'duration',duration,'smooth',smooth,'mode',mode);   
    if mode == 'ccg'
        output(S) = max(ccg(:,1,2)./sum(ccg(:,1,2)));
    elseif mode == 'ccv'
        output(S) = max(ccg(:,1,2));
    end
end
end