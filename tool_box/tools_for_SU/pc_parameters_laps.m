function [within_mean1,within_mean2,between] = pc_parameters_laps(n_lap,pos_1,spks_1,in_lap1,pos_2,spks_2,in_lap2,sigma,Xedges,varargin)
%
% pc_parameters_laps -
% 
% ---> INPUTS
% n_lap: number of laps to subsample group 1
% pos_: matrix, position timestamps in sec
% spks_: vector, spike time in sec
% in_lap: 2 column matrix with beginning and end of each lap 
% sigma: int,FiringCurve parameter for smoothing rate maps
% Xedges: int,FiringCurve parameter to define the nBins of the rate map
%
% ---> OUTPUTS
% 

%other functions:FiringCurve (FMAtoolbox)
%Azul Silva, 2024

mean1=[];
mean2=[];
between = []; 
for c=1:100
    
    % Randomly takes n_laps from group 1
    temp = in_lap1(randperm(size(in_lap1,1)),:); 
    in_lap1 = sortrows(temp(1:n_lap,:),1);
   
    %Keep only spk and pos of the subsampled laps  
    spks_1 = Restrict(spks_1, in_lap1);
    pos_1 = Restrict(pos_1, in_lap1); 
    
    if ~isempty(varargin) && isequal(varargin{1}, 'velocity_norm_ave')
        velocity_norm_ave= varargin{2}; 
    end 
    if ~isempty(varargin) && isequal(varargin{3}, 'velocity_norm_rew')
        velocity_norm_rew= varargin{4}; 
    end 
    
    
    
    %Calculate within pc parameters  
    [within_1] = Within_lap_random(pos_1,spks_1,in_lap1,sigma,Xedges,'velocity_norm',velocity_norm_rew);
    [within_2] = Within_lap_random(pos_2,spks_2,in_lap2,sigma,Xedges,'velocity_norm',velocity_norm_ave);
    
    % Save 
    mean1=[mean1;within_1];
    mean2=[mean1;within_2];
    
    %Calculate between pc parameters
    [between_lap_r] = Between_lap_random(pos_1,spks_1,in_lap1,pos_2,spks_2,in_lap2,sigma,Xedges,'velocity_norm_ave',velocity_norm_ave,'velocity_norm_rew',velocity_norm_rew);
    
    between = [between;between_lap_r]; %save

end 

within_mean1 = nanmean(mean1);
within_mean2 = nanmean(mean2);
between = nanmean(between);


