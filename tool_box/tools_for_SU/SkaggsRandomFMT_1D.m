function q = SkaggsRandomFMT_1D(spks, pos,ini_time, end_time, smooth, Xedges, per)
% This function asign random positions to the spikes 100 times, each time 
% skaggs was calculated and then 90-quantile was defined.
%
% --- Inputs ---
% spks: vector, spikes times
% pos: matrix, positions in x and their respective time stamps (column1: time / column2: position in x axis)

% smooth: int for smooth the firing curve
% Xedges: int for spatial bins for the firing curve
% per: float, define the percentile of interest. (from 0 to 1)
% ini_time: session initial time (msec) 
%end_time: session final time (msec) 
% --- OUTPUT ---
% q: float, 50-quantile from random skaags values.
% Silva Azul 2023 - edited by Morici Juan Facundo

info = zeros(100,1);
for S = 1:100
    
    %op1: shuffle positions, keep spkies 
     spks_R=spks;
     xpos = pos(:,2); 
     xpos= xpos(randperm(length(xpos))); 
     pos = [pos(:,1),xpos]; 
    
%     % op2: create the same # of random time stamps 
%     n_spk = size(spks,1);
%     time_stamps=ini_time:0.0015:end_time; % 1.5msec is the minimum bio time in between spks
%     idx = sort(randi(size(time_stamps,2),[n_spk 1])); 
%     spks_R = time_stamps(idx)'; 
    
    
    [curve , stats] = FiringCurve(pos, spks_R , 'smooth' , smooth , 'nBins' , Xedges , 'minSize' , 6 , 'minPeak' , 0.2);
    info(S,1) = stats.specificity;
end
    q = quantile(info, per);
end 