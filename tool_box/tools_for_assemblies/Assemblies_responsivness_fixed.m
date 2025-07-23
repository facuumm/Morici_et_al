function [curve , bins , responsive] = Assemblies_responsivness_fixed(patterns,cond,SpikeTrain,events,period,window,bin,smooth,normalization,th)
% Rate tuning curve calculation from assemblies using peaks detected above
% 5SD
% This function construct a firing curve sourrounding an event and then, it
% evaluates if the Single-Units (SU) increase or decrease their response.
%
% [curve , bins , responsive] = SU_responsivness(spikes,clusters,events,limits,period,window,bin,smooth,normalization,th)
%
% --- INPUTS ---
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
%
% events: column vector, contains the time stamps of the events
%
% limits: row vector, it contains the begining and end of the period of
%         time used for the tuning curve normalization.
%         (1st column: begining / 2nd column: end)
%
% period: row vector, it contains the begining and end of the period of
%         time used to calculate the mean response of the SU during the event.
%         (1st column: begining / 2nd column: end)
%         If the event started at 0 and finished at +1, then [0 1]
%
% window: float, total duration of the tuning curve.
%
% bin: float, time window (sec) for Spike Train construction
%
% smooth: int, SD for gaussian kernel.
%
% normalization: string, 'zscore', 'gain', 'none'
%                both 'zscore' and 'gain' use the mean firing rate outside 
%                the events. If 'none' or 'gain', it will shuffle the spks
%                100 times and at each iteration it will calculate the mean
%                FR/gain during the period of interest to create a random
%                distribution to use the th and detect responsive cells.
%
% th: float, threshold is SD to define if a SU is responsive or not.
%
% --- OUTPUTS ---
% curve: matrix, it contains the tuning curves for each assemblie.
%        Example:
%                    A1   A2   A3  ...  A10
%                   Re1  Re1  Re1  ...  Re1
%                   Re2  Re2  Re2  ...  Re2
%                   ...  ...  ...  ...  ...
%                   ReX  ReX  ReX  ...  ReX
%
% bins: vector containing the time bins for plotting tuning curves.
%
% responsive: vector containing tags for responsive and non-responsive
%            Assembly.
%             if 1, increased response.
%             if 0, no changment in the response.
%             if -1, decreased response.
%
% Morici Juan Facundo 07/2025
% Other funtions: binspikes, CCG from FMAToolbox


time = SpikeTrain(:,1);
spks = SpikeTrain(:,2:end);

%calculation of assemblies activity
a = assembly_activity(patterns(:,cond) , spks');

%zscored of assemblies activity
a = zscore(a,1,2);

% Generation of structures to store outputs
curve = [];
responsive = [];

for i = 1:size(a,1)
    %detection of peaks
    y = a(i,:)>=2;
    y = time(y);
    x = events;
    
    [s,ids,~] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
    [ccg,tttt] = CCG(s,ids,'binSize',bin,'duration',10,'smooth',1,'mode','ccg'); %ccg calculation
    bins = tttt;
    
    if size(ccg,2)==1
        ccg = nan(length(bins),2,2);
    end
    
    if strcmp(normalization,'zscore')
        final = zscore(((ccg(:,1,2)./bin)./size(events,1)));
    elseif strcmp(normalization,'gain')
        final = ((ccg(:,1,2)./bin)./size(events,1))./m;
    else
        final = (ccg(:,1,2)./bin)./size(events,1);   
    end

    % Responsivness
    if and(strcmp(normalization,'none'), exist('th','var'))
        [h start] = min(abs(bins-period(1)));
        [h stop] = min(abs(bins-period(2)));
        
        if nanmean(final(start:stop)) >= m+st*th
            responsive = [responsive , 1];
        elseif nanmean(final(start:stop)) <= m-st*th
            responsive = [responsive , -1];
        else
            responsive = [responsive , 0];
        end
    end
    
    % Responsivness if none normalization
        if or(strcmp(normalization,'zscore') , strcmp(normalization,'gain'))

            shuffled_ccg = [];
            for ii = 1 : 200 
                yy = ShuffleSpks(y);
                
                [s,ids,groups] = CCGParameters(x,ones(length(x),1),yy,ones(length(yy),1)*2);
                [ccg,tttt] = CCG(s,ids,'binSize',bin,'duration',window,'smooth',0,'mode','ccg'); %ccg calculation
                bins = tttt;
                
                if size(ccg,2)==1
                    ccg = nan(length(bins),2,2);
                end                
                
                %Compute fr and std outside events    
                if strcmp(normalization,'zscore')
                    shuff = zscore(((ccg(:,1,2)./bin)./size(events,1)));
                elseif strcmp(normalization,'gain')
                    shuff = ((ccg(:,1,2)./bin)./size(events,1))./m;
                end
                
                [h start] = min(abs(bins-period(1)));
                [h stop] = min(abs(bins-period(2)));
                shuffled_ccg = [shuffled_ccg; mean(shuff(start:stop))];
                
            end
                    
            th =  quantile(shuffled_ccg,0.9);
            th1  =  quantile(shuffled_ccg,0.1);
            
            if mean(final(start:stop)) >=  th
                responsive = [responsive , 1];
            
            elseif mean(final(start:stop)) <= th1
                responsive = [responsive , -1];
            
            else
                responsive = [responsive , 0];
            end
       end
    
    if smooth > 0
        curve = [curve , Smooth(final,smooth)];
    else
        curve = [curve , final];
    end
    
    clear x y s ids groups ccg tttt is tmp b m st zscored h start stop ist
end

end