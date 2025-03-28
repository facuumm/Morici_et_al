function [curve1] = pc_firigcurve_laps(n_lap,pos_1,spks_1,in_lap1,sigma,Xedges,varargin)
%
%
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

for c=1:100
    
    % Randomly takes n_laps from group 1
    temp = in_lap1(randperm(size(in_lap1,1)),:); 
    in_lap1 = sortrows(temp(1:n_lap,:),1);
   
    %Keep only spk and pos of the subsampled laps  
    spks_1 = Restrict(spks_1, in_lap1);
    pos_1 = Restrict(pos_1, in_lap1); 
   
    [curve1] = FiringCurve(pos_1, spks_1 , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4 , 'minPeak' , 0.2,'minTime',0.15);
   
    if isempty(varargin)
        mean1=[mean1;curve1.rate];
    else % if there is normalization factor 
        if isequal(varargin{1}, 'norm1')
        norm1= varargin{2};mean1=[mean1;curve1.rate./norm1];
        end 
    end 
     
end 

curve1 = nanmean(mean1);


