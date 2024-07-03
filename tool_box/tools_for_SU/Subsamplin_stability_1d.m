% Subsampling + stability 1st and 2nd half


function [sub_stab_ave, sub_stab_rew] = Subsamplin_stability_1d(pos_ave,pos_rew,spks_ave,spks_rew,Xedges,sigma)
    
     sub_stab_ave = [];
     sub_stab_rew = [];
     for xx=1:100
        %%%%%Subsampling%%%%%%
        %Assigne positions to bins for each condition (rew/ave)
        [N_ave,~,binX]  = histcounts(pos_ave(:,2),Xedges);
        [N_rew,~,binXr] = histcounts(pos_rew(:,2),Xedges);

        % Create matrices with bin info 
        pos_ave = [pos_ave,binX];  
        pos_rew = [pos_rew,binXr];

        %Assigne spikes to bins 
        %Create time serie matrix for spikes: 
        bspks_ave = [pos_ave(:,1),zeros(size(pos_ave(:,1),1),1)]; 
        for i= 1:size(bspks_ave,1)-1
          %Define interval time window: 
           ini = bspks_ave(i,1);
           fin = bspks_ave(i+1,1);
          %Select time_spk within that interval
           in= spks_ave(spks_ave>=ini & spks_ave<fin); 
          %Save results 
          bspks_ave(i,2)= length(in);
        end
      
        bspks_rew = [pos_rew(:,1),zeros(size(pos_rew(:,1),1),1)]; 
        for i= 1:size(bspks_rew,1)-1
          %Define interval time window: 
           ini = bspks_rew(i,1);
           fin = bspks_rew(i+1,1);
          %Select time_spk within that interval
           in= spks_rew(spks_rew>=ini & spks_rew<fin); 
          %Save results 
          bspks_rew(i,2)= length(in);
        end

        %Compute differences between occup maps 
        d = N_ave - N_rew;

        %Create subsampled pos ans spk
        for r=1:size(d,2)
          if d(r)>0 % ave higher counts that rew
              ave_bin = find(pos_ave(:,3)==r);
              bin = pos_ave(ave_bin,:); %timestamps in that bin
              k = sort(randperm(numel(ave_bin),N_ave(r)-N_rew(r)));%randomly chose the same amount of timestamps in rew than in ave 
              sub_timestamp = bin(k,1);
              %Subsample rew position
              pos_ave(ismember(pos_ave(:,1),sub_timestamp),:)=[];
              %Subsample rew spk 
              bspks_ave(ismember(bspks_ave(:,1),sub_timestamp),:)=[];
          elseif d(r)<0 %rew higher counts that ave
              rew_bin = find(pos_rew(:,3)==r);
              bin = pos_rew(rew_bin,:); %timestamps in that bin
              k = sort(randperm(numel(rew_bin),N_rew(r)-N_ave(r)));%randomly chose the same amount of timestamps in rew than in ave 
              sub_timestamp = bin(k,1);
              %Subsample rew position
              pos_rew(ismember(pos_rew(:,1),sub_timestamp),:)=[];
              %Subsample rew spk 
              bspks_rew(ismember(bspks_rew(:,1),sub_timestamp),:)=[];
          elseif d(r)==0   % equal counts
                % not subsampling  
          end
                   
        end
   
        sub_pos_ave = pos_ave(:,1:2);
        sub_pos_rew = pos_rew(:,1:2);
      
        sub_spk_ave = bspks_ave(bspks_ave(:,2)~= 0,:);
        sub_spk_ave = repelem(sub_spk_ave(:,1), sub_spk_ave(:,2)); 
 
        sub_spk_rew = bspks_rew(bspks_rew(:,2)~= 0,:);
        sub_spk_rew = repelem(sub_spk_rew(:,1), sub_spk_rew(:,2)); 
        
        %%%%%Stability%%%%%%
        %----Aversive----%
        % 1st half vs 2nd half 
        half_tresh = round((size(sub_pos_ave,1))/2);
             
        % Dividing spks and pos by half_tresh
         pos_1=sub_pos_ave(1:half_tresh,:); 
         pos_2= sub_pos_ave(half_tresh+1:end,:); 
         spks_1 = sub_spk_ave(sub_spk_ave<sub_pos_ave(half_tresh,1));
         spks_2 = sub_spk_ave(sub_spk_ave>sub_pos_ave(half_tresh,1));
             
         %Calculate remapping parameters 
         [curve1,stats1] = FiringCurve(pos_1, spks_1 , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 2 , 'minPeak' , 0.1);
         [curve2,stats2] = FiringCurve(pos_2, spks_2, 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 2 , 'minPeak' , 0.1);
         
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
         if and(sum(curve2.rate)>0,sum(curve1.rate)>0)
            s = corrcoef(curve1.rate, curve2.rate);
            spatial = s(1,2);
         else 
              spatial = nan; 
         end 
    
         %Peak shift 
         shift = abs(stats1.x(1) - stats2.x(1));
            
         sub_stab_ave = [sub_stab_ave;spatial, fr_change, overlap, shift]; 
         
         %----Reward----%
         % 1st half vs 2nd half 
          half_tresh = round((size(sub_pos_rew,1))/2);
             
         % Dividing spks and pos by half_tresh
         pos_1=sub_pos_rew(1:half_tresh,:); 
         pos_2= sub_pos_rew(half_tresh+1:end,:); 
         spks_1 = sub_spk_rew(sub_spk_rew<sub_pos_rew(half_tresh,1));
         spks_2 = sub_spk_rew(sub_spk_rew>sub_pos_rew(half_tresh,1));
             
         %Calculate remapping parameters 
         [curve1,stats1] = FiringCurve(pos_1, spks_1 , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 2 , 'minPeak' , 0.1);
         [curve2,stats2] = FiringCurve(pos_2, spks_2, 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 2 , 'minPeak' , 0.1);
         
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
        if and(sum(curve2.rate)>0,sum(curve1.rate)>0)
            s = corrcoef(curve1.rate, curve2.rate);
            spatial = s(1,2);
         else 
              spatial = nan; 
         end 
    
         %Peak shift 
         shift = abs(stats1.x(1) - stats2.x(1));
            
         sub_stab_rew = [sub_stab_rew;spatial, fr_change, overlap, shift]; 
      
          
      
     end
     sub_stab_ave = nanmean(sub_stab_ave);
     sub_stab_rew = nanmean(sub_stab_rew);
end