function surrogate = surrogate_ccg(times1,times2,win,dur,b,sm,limits)
% This function create a surrogate distribution of maximal value
% within the defined time window (win). It construct the ccg using
% times1 vs shuffled(times2) 200 times
% Facundo Morici, 07/2024
surrogate = [];
if not(isempty(limits))
    times1 = Restrict(times1,limits);
end
for i = 1:500
    shuf = ShuffleSpks2(times2);
%     shuf = Restrict(shuf,limits);
    if and(length(times1)>5 , length(shuf)>5)
        [s,ids,groups] = CCGParameters(times1,ones(length(times1),1),shuf,ones(length(shuf),1)*2);
        [ccg,T] = CCG(s,ids,'binSize',b,'duration',dur,'smooth',sm,'mode','ccg');
        ccg = (ccg(:,1,2)./length(times1))./b;
        [~,up] = min(abs(T-win));    [~,down] = min(abs(T-(-win)));
        surrogate = [surrogate ; max(ccg(down:up))];
    end
end

end