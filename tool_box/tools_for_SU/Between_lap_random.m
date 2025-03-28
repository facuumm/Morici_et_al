function [between_mean] = Between_lap_random(pos_ave,spks_ave,in_lapA,pos_rew,spks_rew,in_lapR,sigma,Xedges,varargin)
%
% Between_pc - Randomly splits laps of the two conditions (ave and rew) in two groups keeping similar #spks in 
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

between = nan(100,4);

for c = 1:100
    %%% Aversive 
    temp = in_lapA(randperm(size(in_lapA,1)),:); 
    half = ceil(size(in_lapA,1)/2); 
   
    g1= sortrows(temp(1:half,:),1);
    g2= sortrows(temp(half+1:end,:),1);
   
    %Split spk and pos from each lap 
    spks_1 = Restrict(spks_ave,g1);
    pos_1 = Restrict(pos_ave,g1); 
    
    spks_2 = Restrict(spks_ave,g2);
    pos_2 = Restrict(pos_ave,g2); 

    ave{1,1}= spks_1; ave{1,2}= pos_1; 
    ave{2,1}= spks_2; ave{2,2}= pos_2; 
    clear  g1 g2 spks_1 pos_1 spks_2 pos_2 temp half 
    
    %%% Rewarded  
    temp = in_lapR(randperm(size(in_lapR,1)),:); 
    half = ceil(size(in_lapR,1)/2); 
   
    g1= sortrows(temp(1:half,:),1);
    g2= sortrows(temp(half+1:end,:),1);
   
    %Split spk and pos from each lap 
    spks_1 = Restrict(spks_rew,g1);
    pos_1 = Restrict(pos_rew,g1); 
    
    spks_2 = Restrict(spks_rew,g2);
    pos_2 = Restrict(pos_rew,g2); 

    rew{1,1}= spks_1; rew{1,2}= pos_1; 
    rew{2,1}= spks_2; rew{2,2}= pos_2; 
    clear  g1 g2 spks_1 pos_1 spks_2 pos_2 temp half 
 
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
    
    if ~isempty(varargin) && isequal(varargin{1}, 'velocity_norm_ave')
        velocity_norm_ave= varargin{2}; 
        curveA.rate= curveA.rate./velocity_norm_ave;
    end 
    if ~isempty(varargin) && isequal(varargin{3}, 'velocity_norm_rew')
        velocity_norm_rew= varargin{4}; 
        curveR.rate= curveR.rate./velocity_norm_rew;
    end 
    
    %Fr change
    fr_change = abs((fr_A - fr_R)/(fr_A + fr_R));

    %Change in fr
    change_fr2 = (fr_R - fr_A)/(fr_R + fr_A);
    
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
    between(c,1)=spatial;
    between(c,2)=fr_change;
    between(c,3)=overlap;
    between(c,4)=shift;
    between(c,5)=change_fr2;

% Control plot 
% figure;hold on; 
% subplot(2,1,1);imagesc(curveA.rate), colormap 'jet'
% subplot(2,1,2);imagesc(curveR.rate), colormap 'jet'
% sgtitle('Firing curvs laps');
%    
end 
between_mean = nanmean(between);
end