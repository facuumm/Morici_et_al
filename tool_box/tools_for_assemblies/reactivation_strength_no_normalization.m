function [R] = reactivation_strength_no_normalization(patterns , cond , SpikeTrain , TS, normalization)
% This function calculate Reactivation Strength (van de Ven et al (2016)).
% Find the peaks of the zscored assemblies strength and calculate the mean
% during Pre sleep and Post sleep. Then, it calculates Post-Pre.
%
% [R] = reactivation_strength(patterns , cond , SpikeTrain , is , th , type)
%
% --- Inputs ---
% patterns: matrix with the weigths for each cell in each assembly.
%           Structure: Single-Units x Assemblies (rows x column)
%
% cond: logical, to select which pattern will be used (1, include / 0, not include)
%       Row vector with the same elements as number of assemblies.
%       Example,
%       patterns:     A1    A2    A3
%                 SU1 0.6   0.1   0.1       ---> cond: 0   1   0
%                 SU2 0.3   0.6   0.1       In this case, only the second
%                 SU3 0.1   0.3   0.3       assembly (A2) will be selected.
%                 SU4 0.1   0.1   0.6
%
% SpikeTrain: matrix, First column, time bins, Rest columns Spike Counts
%             Example:   t    SU1   SU2   SU3   SU4
%                        0     3     5     1     3    This Spike Train was
%                       0.1    1     2     1     0    constructed using a
%                       0.2    3     5     1     3    0.1 sec time bin.
%                       ...   ...   ...   ...   ...
%
% Is: structure, it contains the periods of time (in sec) of interest.
%     Struture:
%                 Is.baseline     Should be a vector with the same length
%                 Is.aversive     of the SpikeTrain with logical values to
%                 Is.reward       include or not the time bin.
%                 Is.runaversive
%                 Is.runreward
%
%     Example:
%             is.aversive:
%                         b1  b2  b3  b4  b5  b6  b7  b8  b9
%                         0   0   1   1   1   1   0   0   0
%
% th: float/int, threshold value for peaks detection. Note that this
%     function zscored the assemblies activity.
%
% type: String, 'A' if is Aversive assemblies 'R' if they are Reward
%
% config: int, if is 1 Aversive occurs first, if is 2, Reward was first
%
% normalization: int, if is 1 then (mean(Post) - mean(Pre))/mean(Post+Pre)
%                     if is 0 then (mean(Post) - mean(Pre))
%
% templates: 2-columns matrix, first column for dHPC and second for vHPC.
%            for more details see assembly_activity_only_joint
%
% --- OUTPUT ---
% R: column vector storing all the Reactivation Strength values for sleep
%    and awake periods.
%
%    S1 ---- *C1*  ---- S2 ----  ---- C2 ---- S3
%
%    * Assemblies detected during this session. It could do the same if C2
%    is the condition of interest. Just be sure you are choosing correctly
%    the type and config inputs.
%
%    Structure:  mean(S1)-mean(S2)   mean(C1)   mean(C2)
%                         R1           A11        A21
%                         R2           A12        A22
%                         R3           A13        A23
%                         ...          ...        ...
%
% requirments:
%       assembly_activity.m from Lopes-dos-Santos et al 2013 (*, see below)
%
% Morici Juan Facundo 02/2024

bins = SpikeTrain(:,1);
dt = bins(2)-bins(1); % delta time
spks = SpikeTrain(:,2:end);
iterations = 1;
prc = 75;

a = assembly_activity(patterns(:,cond) , spks');

A = a;
a = zscore(a,1,2);
R = [];

% Shorten the same amout of time across sleep epochs
Pre = sum(TS.Pre(:,2) - TS.Pre(:,1));
Post = sum(TS.Post(:,2) - TS.Post(:,1));

if Pre < Post
    TS.Post = shortenIntervals(TS.Post,Pre,false);
else
    TS.Pre = shortenIntervals(TS.Pre,Post,true);
end
clear Pre Post

%     TS.Post = shortenIntervals(TS.Post,1200,false);
%     TS.Pre = shortenIntervals(TS.Pre,1200,true);



for i = 1:size(a,1)

    Pre = InIntervals(bins,TS.Pre);
%     Pre = nanmean(a(i,Pre));    
    [Pre,loc] = findpeaks(a(i,Pre),bins(Pre),'MinPeakHeight',5);    
    Pre = nanmean(Pre);
    
    Post = InIntervals(bins,TS.Post);
%     Post = nanmean(a(i,Post));    
    [Post,loc] = findpeaks(a(i,Post),bins(Post),'MinPeakHeight',5);    
    Post = nanmean(Post);   

    
    if normalization
        R = [R ; Post/Pre];
    else
        R = [R ; Pre Post];
    end
    clear Pre Post
end

end