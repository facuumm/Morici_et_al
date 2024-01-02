%Auxiliary function of Between.m and Within.m
%It splits pos(1D) and spks(spk times in sec) into bins of bin_size. Then
%it randomly splits bins into 2 groups of the same size
%Azul Silva, 2023

function [groups] = Bins_half(pos,spks,bin_size)

    %Asigne tspk to position bins 
    tspk = cell(size(pos(:,1),1),1);
    for ind =1:size(spks,1)
        nearest = spks(ind);
        [~,indice]=min(abs(pos(:,1)- nearest));
        if isempty(tspk(indice))
            tspk{indice} = spks(ind);
        else
            tspk{indice} = [tspk{indice},  spks(ind)];                                
        end
    end

    %Calculate the n° of timestamps in bin_size sec 
    sampling_fr = mean(diff(pos(:,1)));
    n_timestamps_bin = ceil(bin_size/sampling_fr);
                    
    %Split pos data and corrsponding tspk in bins
    n = numel(pos(:,2));
    bin_x = mat2cell(pos(:,2),diff([0:n_timestamps_bin:n-1,n]));
    bin_time = mat2cell(pos(:,1),diff([0:n_timestamps_bin:n-1,n]));
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

    %Save output 
    groups{1,1}=time_x_y_1;
    groups{1,2}=tspk_1';
    groups{2,1}=time_x_y_2;
    groups{2,2}=tspk_2';
end