%Subsampling firing rate 1D
%Funcion que binea un time serie y subsamplea para que las dos condiciones
%tengan la misma cantidad de spikes per bin.
%azul
% NOT WORKING 

function [sub_ratemap_ave, sub_ratemap_rew] = Subsampling_fr_1d(pos_ave, spks_ave,pos_rew, spks_rew,dt,Xedges,sigma)

     %Assigne positions to bins for each condition (rew/ave)
     [N_ave,~,~]  = histcounts(pos_ave(:,2),Xedges);
     [N_rew,~,~] = histcounts(pos_rew(:,2),Xedges);

    %Create spk map
     Xs = interp1(pos_ave(:,1),pos_ave(:,2),spks_ave);
     Nspk_ave = histcounts(Xs,Xedges);
     
     Xs = interp1(pos_rew(:,1),pos_rew(:,2),spks_rew);
     Nspk_rew = histcounts(Xs,Xedges);
     
     %Compute differences between spike maps 
     d = Nspk_ave - Nspk_rew;

     %Create a subsampled occup map and Nspike   
      sub_n_ave = nan(size(d));
      sub_n_rew = nan(size(d));
      
      for r=1:size(d,2)
          if d(r)>0 % ave higher counts that rew
               sub_n_rew(r) = Nspk_rew(r);
               sub_n_ave(r) = Nspk_rew(r); %subsampled 
          elseif d(r)<0 %rew higher counts that ave
               sub_n_ave(r) = Nspk_ave(r);
               sub_n_rew(r) = Nspk_ave(r); %subsampled
          elseif d(r)==0   % equal counts
                sub_n_ave(r) = Nspk_ave(r);
                sub_n_rew(r) = Nspk_rew(r);
                % do not touch spk  
          end
                   
      end  
      
      occmap_ave = N_ave*dt; 
      occmap_rew = N_rew*dt; 
      

      %Filtering
      occmap_ave = Smooth(occmap_ave,sigma)'; 
      occmap_rew = Smooth(occmap_rew,sigma)'; 
%       subplot(2,1,2);PlotColorMap(sub_occmap);
      sub_n_ave = Smooth(sub_n_ave,sigma)'; 
      sub_n_rew = Smooth(sub_n_rew,sigma)'; 
                
             
      sub_ratemap_ave = sub_n_ave./occmap_ave;
      sub_ratemap_rew = sub_n_rew./occmap_rew;

end




function [sub_pos_ave,sub_spk_ave,sub_pos_rew,sub_spk_rew] = Subsampling_1d(pos_ave, spks_ave,pos_rew, spks_rew,dt,Xedges,sigma)

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
      
end


