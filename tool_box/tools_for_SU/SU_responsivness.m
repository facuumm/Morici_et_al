function [curve , bins , responsive] = SU_responsivness(spikes,clusters,events,limits,period,window,bin,smooth,normalization,th)
% Fining rate tuning curve calculation.
% This function construct a firing curve sourrounding an event and then, it
% evaluates if the Single-Units (SU) increase or decrease their response.
%
% [curve , bins , responsive] = SU_responsivness(spikes,clusters,events,limits,period,window,bin,smooth,normalization,th)
%
% --- INPUTS ---
% spikes: matrix, spikes times
%         (1st column: cluster id / 2nd column: time stamps)
%
% clusters: column vector, contains the cluster ids
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
% curve: matrix, it contains the tuning curves for each SU.
%        Example:
%                   SU1  SU2  SU3  ...  SU10
%                   Re1  Re1  Re1  ...  Re1
%                   Re2  Re2  Re2  ...  Re2
%                   ...  ...  ...  ...  ...
%                   ReX  ReX  ReX  ...  ReX
%
% bins: vector containing the time bins for plotting tuning curves.
%
% responsive: vector containing tags for responsive and non-responsive SU
%             if 1, increased response.
%             if 0, no changment in the response.
%             if -1, decreased response.
%
% Morici Juan Facundo, 09/2023
% Other funtions: binspikes, CCG from FMAToolbox

curve = [];
responsive = [];
for i = 1 : size(clusters,1)
    x = events;
    y = spikes(spikes(:,1)==clusters(i),2);
    
    [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
    [ccg,tttt] = CCG(s,ids,'binSize',bin,'duration',window,'smooth',0,'mode','ccg'); %ccg calculation
    bins = tttt;
    
    %Compute fr and std outside events 
    ist = InvertIntervals([events-4 events+4],limits(1),limits(2)); % keep intervals outside event
    [tmp,b]=binspikes(y,1/bin,limits);
    tmp = tmp./bin; % binned firing rate
    is = InIntervals(b,ist); % binned firing rate outside events
%     m = length(Restrict(y,ist))/sum(ist(:,2)-ist(:,1)); 
    m = mean(tmp(is),'omitnan'); % mean fr outside events
    st = std(tmp(is),0,1,'omitnan'); % std outside events
    
    if strcmp(normalization,'zscore')
        final = (((ccg(:,1,2)./bin)./size(events,1))-m)./st;
%         figure , plot(bins,final)
            
    elseif strcmp(normalization,'gain')
        final = ((ccg(:,1,2)./bin)./size(events,1))./m;
%         figure , plot(bins,final)
    
    else
        final = (ccg(:,1,2)./bin)./size(events,1);
%         figure , plot(bins,final)
    
    end

    % Responsivness if Zscore or Gain
    if and(strcmp(normalization,'zscore'), exist('th','var'))
        [h start] = min(abs(bins-period(1)));
        [h stop] = min(abs(bins-period(2)));
        
        if mean(final(start:stop)) >= th
            responsive = [responsive , 1];
            
        elseif mean(final(start:stop)) <= th*-1
            responsive = [responsive , -1];
            
        else
            responsive = [responsive , 0];
        end
    end
    
    % Responsivness if none normalization
        if or(strcmp(normalization,'none') , strcmp(normalization,'gain'))
        [h start] = min(abs(bins-period(1)));
        [h stop] = min(abs(bins-period(2)));
        
        threshold = [];
        for ii = 1 : 200 % en este loop estoy corriendo  200 veces lo mismo! 
            x = events;
            yy = ShuffleSpks(y); % ATTENTION YY is not being used! it should be used in line 118
            
            [s,ids,groups] = CCGParameters(x,ones(length(x),1),y,ones(length(y),1)*2);
            [ccg,tttt] = CCG(s,ids,'binSize',bin,'duration',window,'smooth',0,'mode','ccg'); %ccg calculation
            if strcmp(normalization,'none')
                threshold = [threshold ; (ccg(start:stop,1,2)./bin)./size(events,1)];
            elseif strcmp(normalization,'gain')
                threshold = [threshold ; ((ccg(start:stop,1,2)./bin)./size(events,1))./m];
            end
        end
        
        if mean(final(start:stop)) >= std(threshold)*th
            responsive = [responsive , 1];
            
        elseif mean(final(start:stop)) <= std(threshold)*th*-1
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