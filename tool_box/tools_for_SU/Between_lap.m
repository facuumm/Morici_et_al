function [between_lap] = Between_lap(pos_ave,spks_ave,in_lapA,pos_rew,spks_rew,in_lapR,sigma,Xedges)
%
% Between_pc - Splits laps of the two conditions (ave and rew) in two groups keeping similar #spks in 
% each group and randomly selects one group of each condition and calcualtes remapping parameters
%(spatial correlation, firing rate change,rate overlap and peak shift)between conditions. 
%
% [between] = Between_lap(pos_ave,spks_ave,in_lapA,pos_rew,spks_rew,in_lapR,sigma,Xedges)
%
% ---> INPUTS
% pos: matrix, timestamps in sec and position, one variable for each condition (ave and rew)
% spks: vector, spike time in sec, one variable for each condition (ave and rew)
% in_lap: 2 column matrix with beginning and end of each lap (ave and rew)
% sigma: int,FiringCurve parameter for smoothing rate maps
% Xedges: int,FiringCurve parameter to define the nBins of the rate map
%
% ---> OUTPUTS
% between: vector with the remapping parameters(spatial correlation, firing rate change,
%          rate overlap and peak shift).
% 
%
%other functions:FiringCurve (FMAtoolbox), Bins_half
%
%Azul Silva, 2024

%%% Aversive 
% Split laps keeping similar #spks per lap  
%Compute the # spk per lap
[~,interval,~] = InIntervals(spks_ave,in_lapA);
    
lap_spk = [(1:size(in_lapA,1))',in_lapA]; 
for i=1:size(lap_spk,1)
    lap_spk(i,4) = sum(interval==i); 
end
clear i    

%Split laps 
g =1; g1= []; g2=[];
while ~isempty(lap_spk)
    [~,idx] = max(lap_spk(:,4)); 
    if g==1
        g1=[g1;lap_spk(idx,:)];
        g = 2;
    elseif g==2
        g2=[g2;lap_spk(idx,:)];
        g = 1;   
    end
    lap_spk(idx,:)=[]; 
end
clear lap_spk g 

g1= sortrows(g1,1); 
g2= sortrows(g2,1); 

%Split spk and pos from each group
spks_1 = Restrict(spks_ave,g1(:,2:3));
pos_1 = Restrict(pos_ave,g1(:,2:3)); 
    
spks_2 = Restrict(spks_ave,g2(:,2:3));
pos_2 = Restrict(pos_ave,g2(:,2:3)); 
clear g1 g2    

ave{1,1}= spks_1; ave{1,2}= pos_1; 
ave{2,1}= spks_2; ave{2,2}= pos_2; 
clear  g1 g2 spks_1 pos_1 spks_2 pos_2
    
%%% Rewarded  
% Split laps keeping similar #spks per lap  
%Compute the # spk per lap
[~,interval,~] = InIntervals(spks_rew,in_lapR);
    
lap_spk = [(1:size(in_lapR,1))',in_lapR]; 
for i=1:size(lap_spk,1)
    lap_spk(i,4) = sum(interval==i); 
end
clear i    

%Split laps 
g =1; g1= []; g2=[];
while ~isempty(lap_spk)
    [~,idx] = max(lap_spk(:,4)); 
    if g==1
        g1=[g1;lap_spk(idx,:)];
        g = 2;
    elseif g==2
        g2=[g2;lap_spk(idx,:)];
        g = 1;   
    end
    lap_spk(idx,:)=[]; 
end
clear lap_spk g 

%Split spk and pos from each group
spks_1 = Restrict(spks_rew,g1(:,2:3));
pos_1 = Restrict(pos_rew,g1(:,2:3)); 
    
spks_2 = Restrict(spks_rew,g2(:,2:3));
pos_2 = Restrict(pos_rew,g2(:,2:3)); 
clear g1 g2    

rew{1,1}= spks_1; rew{1,2}= pos_1; 
rew{2,1}= spks_2; rew{2,2}= pos_2; 
clear  g1 g2 spks_1 pos_1 spks_2 pos_2
 
%Choose randomly one of the groups of ave and rew  for the comparison
id = round(rand)*2; 
if id==0 % not the best way
   id =1;
end 
    
time_x_ave = ave{id,2};
tspk_ave = ave{id,1};
    
time_x_rew = rew{id,2};
tspk_rew = rew{id,1};
    
%Calculate remapping parameters 
[curveA, statsA] = FiringCurve(time_x_ave, tspk_ave , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4 , 'minPeak' , 0.2,'minTime',0.15);
[curveR, statsR] = FiringCurve(time_x_rew, tspk_rew , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4 , 'minPeak' , 0.2,'minTime',0.15);
                    
fr_A= nanmean(curveA.rate);
fr_R= nanmean(curveR.rate);
                    
%Fr change
fr_change = abs((fr_A - fr_R)/(fr_A + fr_R));

%Rate overlap
if fr_A<=fr_R 
    overlap = fr_A/fr_R;
else 
    overlap = fr_R/fr_A;
end
    
%Spatial  corr
s = corrcoef(curveA.rate, curveR.rate);
spatial = s(1,2);
    
% Peak shift 
shift = abs(statsA.x(1) - statsR.x(1)); 

%Save 
between_lap = [spatial, fr_change, overlap,shift];    

% Control plot 
% figure;hold on; 
% subplot(2,1,1);imagesc(curveA.rate), colormap 'jet'
% subplot(2,1,2);imagesc(curveR.rate), colormap 'jet'
% sgtitle('Firing curvs laps');
%    

end