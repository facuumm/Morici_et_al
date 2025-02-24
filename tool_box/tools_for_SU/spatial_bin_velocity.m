%Function that spatially bin instantaneous velocity
%INPUT
% data   matrix with timestamps, position, instantaneous velocity 
%Xedges Number of spatial bins
%OUTPUT
%spatial_binned_velocity mean instantaneous velocity per spatial bin

function spatial_binned_velocity =  spatial_bin_velocity(data, Xedges)
    
       %Assigne positions to bins for each condition 
        [~,~,binX]  = histcounts(data(:,2),Xedges);
       % Create matrices with bin info 
      data = [data,binX];  
       %Computed the mean velocity per each spatial bin Xedges
       spatial_binned_velocity = zeros(1, Xedges); 
        for i = 1:Xedges
            temp_velocity = data(data(:,4) == i, 3); % Filter velocitys of bin i
            if ~isempty(temp_velocity) % check if there is data
                 spatial_binned_velocity(i) = nanmean(temp_velocity); % compute the mean 
            else
               spatial_binned_velocity(i) = NaN;
            end
        end
end