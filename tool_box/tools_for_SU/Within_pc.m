function [within_mean,within_percentil ] = Within_pc(pos_tmp,spks_tmp,bin_size,sigma,Xedges)
%
% Within_pc - Randomly splits spike and position data in two and calcualtes remapping parameters
% (spatial correlation, firing rate change, rate overlap and peak shift) within them. Repeat 500 times  
%
% ---> INPUTS
% pos_tmp: matrix, position timestamps in sec
% spks_tmp: vector, spike time in sec
% bin_size: bin size in sec (1 recommended). Time window to binnarized pos
%           and spike data
% sigma: int,FiringCurve parameter for smoothing rate maps
% Xedges: int,FiringCurve parameter to define the nBins of the rate map
%
% ---> OUTPUTS
% within_mean: vector with the mean value of the 500 iterations for the 4
%              remapping parameters
% within_percentil: vector with the percentil value of the 500 iterations for the 4
%                   remapping parameters - to caluclate remappings tresholds;
%
%other functions:FiringCurve (FMAtoolbox)
%
%Azul Silva, 2023
within = nan(500,4);

for c=1:500

    %Asigne tspk to position bins 
    tspk = cell(size(pos_tmp(:,1),1),1);
    for ind =1:size(spks_tmp,1)
        nearest = spks_tmp(ind);
        [~,indice]=min(abs(pos_tmp(:,1)- nearest));
        if isempty(tspk(indice))
            tspk{indice} = spks_tmp(ind);
        else
            tspk{indice} = [tspk{indice},  spks_tmp(ind)];                                
        end
    end

    %Calculate the n° of timestamps in bin_size sec 
    sampling_fr = mean(diff(pos_tmp(:,1)));
    n_timestamps_bin = ceil(bin_size/sampling_fr);
                    
    %Split pos data and corrsponding tspk in bins
    n = numel(pos_tmp(:,2));
    bin_x = mat2cell(pos_tmp(:,2),diff([0:n_timestamps_bin:n-1,n]));
    bin_time = mat2cell(pos_tmp(:,1),diff([0:n_timestamps_bin:n-1,n]));
    bin_spk = mat2cell(tspk,diff([0:n_timestamps_bin:n-1,n]));

    %Split bins randomly into two equal size groups 

    numBins = size(bin_x,1);
    halfNumBins = ceil(numBins/2); %size of each group

    randomIndex = randperm(numBins); %random index to split 

    %Indexs for each group               
    group1Index = randomIndex(1:halfNumBins);
    group2Index = randomIndex(halfNumBins+1:end);  

    %Asigne bins to groups 
    g1_x= bin_x(group1Index);
    g2_x= bin_x(group2Index);

                     
    g1_time= bin_time(group1Index);
    g2_time= bin_time(group2Index);
                    
    g1_spike= bin_spk(group1Index);
    g2_spike= bin_spk(group2Index);

    % Convert cell to vector 
    x_pos_1 = cell2mat(g1_x);
    x_pos_2 = cell2mat(g2_x);
                   
                  
    time_1 = cell2mat(g1_time);
    time_2 = cell2mat(g2_time);

    tspk_1 = sort(flattenCellArray(g1_spike));
    tspk_2 = sort(flattenCellArray(g2_spike));
             
    %Save in matrix 
    time_x_y_1= sortrows([time_1,x_pos_1],1);
    time_x_y_2= sortrows([time_2,x_pos_2],1);

    %Calculate remapping parameters 
    [curve1,stats1] = FiringCurve(time_x_y_1, tspk_1' , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4 , 'minPeak' , 0.2);
    [curve2,stats2] = FiringCurve(time_x_y_2, tspk_2' , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4 , 'minPeak' , 0.2);
                    
    fr_1= nanmean(curve1.rate);
    fr_2= nanmean(curve2.rate);
                    
    %Fr change
    fr_change = abs((fr_1 - fr_2)/(fr_1 + fr_2));

    %Rate overlap
    if fr_1<=fr_2 
        overlap = fr_1/fr_2;
    else 
        overlap = fr_2/fr_1;
    end
    
    %Spatial  corr
    s = corrcoef(curve1.rate, curve2.rate);
    spatial = s(1,2);
    
    %Peak shift 
    shift = abs(stats1.x(1) - stats2.x(1)); 
    
    %Save 
    within(c,1)=spatial;
    within(c,2)=fr_change;
    within(c,3)=overlap;
    within(c,4)=shift;
    
   
end 

within_mean = nanmean(within);

spatial = quantile(within(:,1), 0.1);
fr = quantile(within(:,2), 0.9);
overlap = quantile(within(:,3), 0.1);
shift = quantile(within(:,4), 0.9);
within_percentil = [spatial, fr, overlap,shift]; 

end