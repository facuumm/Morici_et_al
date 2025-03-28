function result = different_from_surrogate(times1,times2,win,dur,b,sm,limits)
% This function create a surrogate distribution of maximal value
% within the defined time window (win). It construct the ccg using
% times1 vs shuffled(times2) 200 times.
% Then it gives you True or False deppending if the maximal value coming
% from the real CCG is different from the 99th percentile of the surrogate.
% 
% --- Output ---
% result:
% [True False] if times 2 happend frequently before times1
% [False True] if times 2 happend frequently after times1
% [False False] none of the previous conditions
%
% Facundo Morici, 03/2025

surrogate1 = [];
surrogate2 = [];

if not(isempty(limits))
    times1 = Restrict(times1,limits);
end

[s,ids,groups] = CCGParameters(times1,ones(length(times1),1),times2,ones(length(times2),1)*2);
[ccg,T] = CCG(s,ids,'binSize',b,'duration',dur,'smooth',sm,'mode','ccg');
ccg = (ccg(:,1,2)./length(times1))./b;

% Definition of time windows to check
[~,up] = min(abs(T-win));    [~,down] = min(abs(T-(-win)));    [~,mid] = min(abs(T-0));

% value1 = max(ccg(down : mid));
value2 = max(ccg(mid : up));


for i = 1:200
    shuf = ShuffleSpks2(times2);
    if and(length(times1)>5 , length(shuf)>5)
        [s,ids,groups] = CCGParameters(times1,ones(length(times1),1),shuf,ones(length(shuf),1)*2);
        [ccg,T] = CCG(s,ids,'binSize',b,'duration',dur,'smooth',sm,'mode','ccg');
        ccg = (ccg(:,1,2)./length(times1))./b;

        surrogate1 = [surrogate1 ; max(ccg(down : mid))];
        surrogate2 = [surrogate2 ; max(ccg(mid : up))];
    end
end

% result = [value1 > prctile(surrogate1,99) value2 > prctile(surrogate2,99)];
result = [value2 > prctile(surrogate2,99)];
end