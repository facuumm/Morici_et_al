function [Pr, prMax] = bayesian_replay(Spks , ids , RateMap , time , dt , d)
% Replay detection and significance testing.
% This function calculates likelihood R between the replayed trajectory and
% a constant-velocity trajectory inferred with the initial and final
% positions. To determine if R is higher than expected by chance, this
% function constructs two distributions of Rs values by:
%       1) Circularly shifting the estimate position at each time bin by a
%          random distance 1000 times.
%       2) Randomly permuting tuning curves 1000 times.
%
% Method developed by Davidson et al 2009, doi:10.1016/j.neuron.2009.07.027
%
% USAGE
%       R = bayesian_replay(RateMap, nSpks, start, stop)
%
% --- INPUTS ---
% Spks: Matrix, first column must contain cluster id, second column times
%
% ids: column vector, cluster ids
%
% RateMap: Matrix, rate map for all the cells across the space
%          (rows: cell ids / columns: spatial bins)
%
% nSpks: column vector, number of spks ocurring during the event for each
%        unit in the RateMap. ** length(nSpks) == length(RateMap) **
%
% time: row vector, contain the start and end of the putative-replay event
%
% dt: float, time interval for bining.
%
% d: float, Distance margin of acceptance
%
% --- OUTPUTS ---
% Pr = [nTemporalBin X nSpatialBins] matrix of posterior probabilities
% prMax = the spatial bin with higher spatial probabilities for each
%       temporalBin in Cr (note, ties go to the lower numbered bins as 
%       consistent with the behavior of the second output of built in function
%       'max')
% Morici Juan Facundo, 09/06/2023
% 116-46
% from Andres Grosmark 2015


% Variables definition
time = [time(1):dt:time(2)];
d = d/3;

% Count Spks per time
Cr = [];
for i = 1:size(time,2)-1
    p = count_spks(Spks, ids , time(i), time(i+1));
    Cr = [Cr ; p'];
end
Cr = Cr.*dt; clear i

[Pr, prMax] = placeBayes(Cr, RateMap, dt); 
% [Pr, prMax] = placeBayesLogBuffered1(Cr, RateMap, dt); 
 
end
% 
% tmp = radon(Pr);
% 
% 
% 
% [i] = max(tmp,[],'all')
% x1 = [1:60];
% x1 = x1(logical(sum(tmp==i,2)));
% 
% V = [1:size(tmp,2)];
% V = V(logical(sum(tmp==i,1)));
% % imagesc(time,[1:60],Pr)
% 
% [r x1] = max(Pr(:,1)); %initial position
% [r x2] = max(Pr(:,end)); clear r  % final position
% % V = (x2 - x1) / (time(end) - time(2)); % constant velocity calculation
% V = [-60 : 1 : 60];
% 
% estimated = [];
% % Construction of estimated trayectory using V
% for i = 1:size(V,2)
%     v = V(i);
%     ee = [];
%     for ii = 1:60
%         e = [];
%         for iii = 1 : size(time,2)-1
%             if iii == 1
%                 e = [e , ii + (v*0*dt)];
%             else
%                 e = [e , ii + (v*iii*dt)];
%             end
%         end
%         ee = [ee ; e];
%     end
%     estimated(:,:,i) = ee;
%     clear v ee
% end
% 
% % template constructed detecting max(Pr) in time
% [r , template] = max(Pr);
% 
% % Replay score calculation
% R = [];
% for i = 1:size(V,2)
%     for ii = 1:60
%         t = (sum(abs(template - estimated(ii,:,i))<=d)/size(template,2))*100;
%     end
%     R(ii,i) = t;
%     clear t
% end
% scatter([1:121] , R(60,:))
% % R = (sum(abs(template - estimated)<=d)/size(template,2))*100; 
% clear template r
% 
% %% Surrogated distributions
% % Circularly shifting the estimate position
% distribution1 = [];
% for i = 1 : 1000
%     r = round(1 + (60-1).*rand(6,1))';
%     tmp = circshift(Pr,r);
%     
%     [r x1] = max(tmp(:,1)); %initial position
%     [r x2] = max(tmp(:,end)); clear r  % final position
%     V = (x2 - x1) / (time(end) - time(1)); % constant velocity calculation
%     
%     % Construction of estimated trayectory using V
%     estimated = [];
%     for i = 1 : size(time,2)-1
%         if i == 1
%             estimated = [estimated , x1 + (V*0*dt)];
%         else
%             estimated = [estimated , x1 + (V*i*dt)];
%         end
%     end
%     
%     % template constructed detecting max(Pr) in time
%     [r , template] = max(tmp);
%     
%     % Replay score calculation
%     distribution1 = [distribution1 ; (sum(abs(template - estimated)<=d)/size(template,2))*100];
%     
% end
% [h1 p1] = lillietest(distribution1,'MCTol',0.01);
% 
% 
% % Randomly permuting tuning curves
% distribution2 = [];
% for i = 1 : 1000
%     tmp = [];
%     xr = RateMap(randperm(size(RateMap,1)),:);
%     for i = 1:size(time,2)-1
%         p = count_spks(Spks, ids , time(i), time(i+1));
%         t = time(i+1) - time(i);
%         f = xr .^ p;
%         y = prod(f , 1);
%         y = y .* exp(-t*sum(xr));
%         c = sum(y);
%         tmp = [tmp , (y ./ c)'];
%         clear p t f y c
%     end
%     
%     [r x1] = max(tmp(:,1)); %initial position
%     [r x2] = max(tmp(:,end)); clear r  % final position
%     V = (x2 - x1) / (time(end) - time(1)); % constant velocity calculation
%     
%     % Construction of estimated trayectory using V
%     estimated = [];
%     for i = 1 : size(time,2)-1
%         if i == 1
%             estimated = [estimated , x1 + (V*0*dt)];
%         else
%             estimated = [estimated , x1 + (V*i*dt)];
%         end
%     end
%     
%     % template constructed detecting max(Pr) in time
%     [r , template] = max(tmp);
%     
%     % Replay score calculation
%     distribution2 = [distribution2 ; (sum(abs(template - estimated)<=d)/size(template,2))*100];
% end
% [h2 p2] = lillietest(distribution2,'MCTol',0.01);
% 
% 
% end