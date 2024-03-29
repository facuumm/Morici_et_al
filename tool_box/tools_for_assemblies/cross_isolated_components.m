function [cross lags] = cross_isolated_components(patterns , cond , SpikeTrain , Is , th , type , config , templates)
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

% if isfield(templates,'dHPC')
%     a = assembly_activity_only_joint(patterns(:,cond) , spks',templates.dHPC,templates.dHPC);
% elseif isfield(templates,'vHPC')
%     a = assembly_activity_only_joint(patterns(:,cond) , spks',templates.vHPC,templates.vHPC);
% elseif isfield(templates,'dvHPC')
a1 = assembly_activity_only_joint(patterns(:,cond) , spks',templates(:,1),templates(:,1));
a2 = assembly_activity_only_joint(patterns(:,cond) , spks',templates(:,2),templates(:,2));
% end

a1 = zscore(a1,1,2);
a2 = zscore(a2,1,2);

cross.pre = [];
cross.post = [];
lags = [];
for i = 1:size(a1,1)
    
%     % Detectioon of peaks
%     [pks1.baseline,loc1.baseline] = findpeaks(a1(i,Is.baseline),bins(Is.baseline));
%     [pks1.reward,loc1.reward] = findpeaks(a1(i,Is.reward),bins(Is.reward));
%     [pks1.aversive,loc1.aversive] = findpeaks(a1(i,Is.aversive),bins(Is.aversive));
%     
%     [pks2.baseline,loc2.baseline] = findpeaks(a2(i,Is.baseline),bins(Is.baseline));
%     [pks2.reward,loc2.reward] = findpeaks(a2(i,Is.reward),bins(Is.reward));
%     [pks2.aversive,loc2.aversive] = findpeaks(a2(i,Is.aversive),bins(Is.aversive));
%     
%     x = loc1.baseline;
%     y = loc2.baseline;
%     [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
%     [ccg,tttt] = CCG(s,ids,'binSize',0.005,'duration',0.4,'smooth',1,'mode','ccg');
%     CB = ccg(:,1,2)./sum(ccg(:,1,2)); clear s ids groups ccg tttt
%     
%     x = loc1.aversive;
%     y = loc2.aversive;
%     [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
%     [ccg,tttt] = CCG(s,ids,'binSize',0.005,'duration',0.4,'smooth',1,'mode','ccg');
%     CA = ccg(:,1,2)./sum(ccg(:,1,2)); clear s ids groups ccg tttt
%     
%     x = loc1.reward;
%     y = loc2.reward;
%     [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
%     [ccg,tttt] = CCG(s,ids,'binSize',0.005,'duration',0.4,'smooth',1,'mode','ccg');
%     CR = ccg(:,1,2)./sum(ccg(:,1,2)); clear s ids groups ccg
%     
    % Restriction of data to events
    a1B =  a1(i,Is.baseline);   a2B =  a2(i,Is.baseline);
    a1R =  a1(i,Is.reward);     a2R =  a2(i,Is.reward);
    a1A =  a1(i,Is.aversive);   a2A =  a2(i,Is.aversive);
    
    % cross-Corr during events
    % during sleep
    [CB lags] = xcorr(a1B , a2B , 20 , 'normalized');
    [CA lags] = xcorr(a1A , a2A , 20 , 'normalized');
    [CR lags] = xcorr(a1R , a2R , 20 , 'normalized');
    
    %     [pks.runreward,loc.reward] = findpeaks(a(i,Is.runreward),bins(Is.runreward),'MinPeakHeight',th);
    %     [pks.runaversive,loc.aversive] = findpeaks(a(i,Is.runaversive),bins(Is.runaversive),'MinPeakHeight',th);
    
    %% Calculation of Reactivation
    if type == 'A' % check if is aversive assembly
        if config == 1
            cross.pre = [cross.pre ; CB];
            cross.post = [cross.post ; CA];
            lags = lags*dt;
%             lags = tttt; clear tttt
        else
            cross.pre = [cross.pre ; CR];
            cross.post = [cross.post ; CA];
            lags = lags*dt;
%             lags = tttt; clear tttt
        end
        
    elseif type == 'R' % check if is reward assembly
        if config == 2
            cross.pre = [cross.pre ; CB];
            cross.post = [cross.post ; CR];
            lags = lags*dt;
%             lags = tttt; clear tttt
        else
            cross.pre = [cross.pre ; CA];
            cross.post = [cross.post ; CR];
            lags = lags*dt;
%             lags = tttt; clear tttt
        end
    end
end

end