function q = SkaggsRandomFMT_1D(spks, pos, smooth, Xedges, per)
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
    
    %op1: circular permutation of position
    spks_R=spks;
    
    trial_duration = pos(end,1)-pos(1,1);
    random_shift =1 + (trial_duration - 1) * rand(); %from 1 sec to the lenght of the session
    %Convert shift to index units
    dt = mean(diff(pos(:,1))); 
    shift_idx = round(random_shift / dt);
    
    permuted_pos= circshift(pos(:,2), shift_idx);
    pos = [pos(:,1),permuted_pos]; 
    
    
    
%     %op3: shuffle positions, keep spkies 
%      spks_R=spks;
%      xpos = pos(:,2); 
%      xpos= xpos(randperm(length(xpos))); 
%      pos = [pos(:,1),xpos]; 
    
%     % op3: create the same # of random time stamps 
%     n_spk = size(spks,1);
%     time_stamps=ini_time:0.0015:end_time; % 1.5msec is the minimum bio time in between spks
%     idx = sort(randi(size(time_stamps,2),[n_spk 1])); 
%     spks_R = time_stamps(idx)'; 
    
    
    [curve , stats] = FiringCurve(pos, spks_R , 'smooth' , smooth , 'nBins' , Xedges , 'minSize' , 6 , 'minPeak' , 0.2);
    info(S,1) = stats.specificity;
end
    q = quantile(info, per);
end 