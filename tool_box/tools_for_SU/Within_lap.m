function [within_lap] = Within_lap(pos_tmp,spks_tmp,in_lap,sigma,Xedges)
%
% Within_lap - Splits laps in two keeping similar # spk in each group. Calcualtes remapping parameters
% (spatial correlation, firing rate change, rate overlap and peak shift) within the two groups.  
%
% ---> INPUTS
% pos_tmp: matrix, position timestamps in sec
% spks_tmp: vector, spike time in sec
% in_lap: 2 column matrix with beginning and end of each lap 
% sigma: int,FiringCurve parameter for smoothing rate maps
% Xedges: int,FiringCurve parameter to define the nBins of the rate map
%
% ---> OUTPUTS
% within_lap: vector with spatial correlation, fr_change, overlap, shift 
%
%other functions:FiringCurve (FMAtoolbox)
%Azul Silva, 2024
  
%% Split laps keeping similar #spks per lap  
%Compute the # spk per lap
[~,interval,~] = InIntervals(spks_tmp,in_lap);
    
lap_spk = [(1:size(in_lap,1))',in_lap]; 
for i=1:size(lap_spk,1)
    lap_spk(i,4) = sum(interval==i); 
end
clear i    
% Split laps 
g =1; g1= []; g2=[];
while ~isempty(lap_spk)
    [~,idx] = max(lap_spk(:,4)); 
    if g==1
        g1=[g1;lap_spk(idx,:)];
        g = 2;
    elseif g==2
        g2=[g2;lap_spk(idx,:)];
        g = 1;   
    end
    lap_spk(idx,:)=[]; 
end

g1= sortrows(g1,1); 
g2= sortrows(g2,1); 

clear lap_spk g 
%% Split spk and pos from each group
spks_1 = Restrict(spks_tmp,g1(:,2:3));
pos_1 = Restrict(pos_tmp,g1(:,2:3)); 
    
spks_2 = Restrict(spks_tmp,g2(:,2:3));
pos_2 = Restrict(pos_tmp,g2(:,2:3)); 
clear g1 g2    
%% Calculate remapping parameters 
[curve1,stats1] = FiringCurve(pos_1, spks_1 , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4 , 'minPeak' , 0.2,'minTime',0.15);
[curve2,stats2] = FiringCurve(pos_2, spks_2 , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4 , 'minPeak' , 0.2,'minTime',0.15);
                    
fr_1= nanmean(curve1.rate);
fr_2= nanmean(curve2.rate);
                    
%Fr change
fr_change = abs((fr_1 - fr_2)/(fr_1 + fr_2));

change_fr2 = (fr_1 - fr_2)/(fr_1 + fr_2);

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
    
% Control plot 
% figure;hold on; 
% subplot(2,1,1);imagesc(curve1.rate), colormap 'jet'
% subplot(2,1,2);imagesc(curve2.rate), colormap 'jet'
% sgtitle('Firing curvs laps');
    

within_lap = [spatial, fr_change, overlap,shift,change_fr2]; 

end




% Extra option: split laps randomly 
 %Randomly splits laps 
%     temp = in_lap(randperm(size(in_lap,1)),:); 
%     half = ceil(size(in_lap,1)/2); 
%    
%     g1= temp(1:half,:);
%     g2=temp(half+1:end,:);
%    
%     %Split spk and pos from each lap 
%     spks_1 = Restrict(spks_tmp,g1);
%     pos_1 = Restrict(pos_tmp,g1); 
%     
%     spks_2 = Restrict(spks_tmp,g2);
%     pos_2 = Restrict(pos_tmp,g2); 
    