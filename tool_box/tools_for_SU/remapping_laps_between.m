function [between_corr] = remapping_laps_between(in_lap_1,pos_1,spks_1,in_lap_2, pos_2,spks_2, sigma, Xedges,min_size, min_peak, min_time)

% Computes spatial correlation across laps of two conditions. Subsamples laps of the condition with more laps. 
% ---> INPUTS
%in_lap_: 2 column matrix with beginning and end of each lap. Spks and pos will be restricted to these laps 
% pos_: matrix, position timestamps in sec
% spks_: vector, spike time in sec
% sigma: int,FiringCurve parameter for smoothing rate maps
% Xedges: int,FiringCurve parameter to define the nBins of the rate map
%
%other functions:FiringCurve (FMAtoolbox)
%Azul Silva, 2024
between_corr = [];
for c=1:100
    %%% Subsampling %%%
    if size(in_lap_1,1)>size(in_lap_2,1)
        n_lap=size(in_lap_2,1);
        % Randomly takes n_laps2 from group 1
        temp = in_lap_1(randperm(size(in_lap_1,1)),:); 
        in_lap1 = sortrows(temp(1:n_lap,:),1);
        in_lap2=in_lap_2;
    elseif size(in_lap_2,1)>size(in_lap_1,1) 
        n_lap=size(in_lap_1,1);
        % Randomly takes n_laps1 from group 2
        temp = in_lap_2(randperm(size(in_lap_2,1)),:); 
        in_lap2 = sortrows(temp(1:n_lap,:),1);
        in_lap1=in_lap_1;
    else 
        in_lap1=in_lap_1;
        in_lap2=in_lap_2;
    end
    
    %%% Computing rate maps per lap %%%
    laps_map1 = [];
    for l=1:size(in_lap1,1)
        %Keep only spk and pos from one
         spks = Restrict(spks_1, in_lap1(l,:));
         pos = Restrict(pos_1, in_lap1(l,:)); 
         [curve1] = FiringCurve(pos, spks , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , min_size , 'minPeak' , min_peak,'minTime',min_time);
         laps_map1 = [laps_map1; curve1.rate];% saving stack of lap rate maps 
    end
    
    laps_map2 = [];
    for l=1:size(in_lap2,1)
        %Keep only spk and pos from one
         spks = Restrict(spks_2, in_lap2(l,:));
         pos = Restrict(pos_2, in_lap2(l,:)); 
         [curve2] = FiringCurve(pos, spks , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , min_size , 'minPeak' , min_peak,'minTime',min_time);
         laps_map2 = [laps_map2; curve2.rate];% saving stack of lap rate maps 
    end
    
%     figure;subplot(1,2,1);imagesc([0:20:140], [1:1:size(laps_map1,1)],laps_map1), colormap 'jet'; title('1');
%     subplot(1,2,2);imagesc([0:20:140], [1:1:size(laps_map2,1)],laps_map2), colormap 'jet'; title('2');
    %%% Computing between spatial corr %%%
    between_corr_matrix = nan(10, 10);
    for i = 1:size(in_lap1,1)
        for j = 1:size(in_lap1,1)
            % Correlate the i-th row of A with the j-th row of B
            between_corr_matrix(i, j) = corr(laps_map1(i, :)', laps_map2(j, :)');
        end
    end
    
    between_corr(:,:,c)=between_corr_matrix;
end    

end 