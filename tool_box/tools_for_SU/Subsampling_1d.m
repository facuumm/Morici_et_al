%Subsampling 1D
%Funcion que binea un time serie y subsamplea para que las dos condiciones
%tengan el mismo muestreo.
%azul
%

function [pos_sub, spikes_sub] = Subsampling_1d(pos_ave, spks_ave,pos_rew, spks_rew,dt, Xedges)

     %Assigne positions to bins for each condition (rew/ave)
     [N_ave,Xedges,binX]  = histcounts(pos_ave(:,2),Xedges);
     [N_rew,Xedgesr,binXr] = histcounts(pos_rew(:,2),Xedges);

     % Create matrices with bin info 
      pos_ave = [pos_ave,binX];  
      pos_rew = [pos_rew,binXr];

     %Assigne spikes to bins 
     %Create time serie matries for spikes: 
      bspks_ave = [pos_ave(:,1),zeros(size(pos_ave(:,1),1),1)]; 
      for i= 1:size(spks_ave,1)-1
          %Define interval time window: 
           ini = bspks_ave(i,1);
           fin = bspks_ave(i+1,1);
          %Select time_spk within that interval
           in= spks_ave(spks_ave>=ini & spks_ave<fin); 
          %Save results 
          bspks_ave(i,2)= length(in);
      end
      
      bspks_rew = [pos_rew(:,1),zeros(size(pos_rew(:,1),1),1)]; 
      for i= 1:size(spks_rew,1)-1
          %Define interval time window: 
           ini = bspks_rew(i,1);
           fin = bspks_rew(i+1,1);
          %Select time_spk within that interval
           in= spks_rew(spks_rew>=ini & spks_rew<fin); 
          %Save results 
          bspks_rew(i,2)= length(in);
     end
              

      %Compute differences between occup maps 
%       d = N_tr - N_ts;
      d = N_ave - N_rew;

     %Create a subsampled occup map and t spike time serie  
      sub_occ = nan(size(d)); 
      for r=1:size(d,1)
          if d(r)>0 % ave higher counts that rew
               sub_occ(r)=N_rew(r); %keep rew counts 
               %Subsampling ave spk 
               %randomly chose the same amount of timestamps in that bin in ave than in rew 
                ave_bin = find(pos_ave(:,3)==r);
                bin = pos_ave(ave_bin,:);
                k = sort(randperm(numel(ave_bin),N_ave(r)-N_rew(r)));
                sub_timestamp = bin(k,1); 
                %Remove those bin time stamp from  bspks_ave
                bspks_ave(ismember(bspks_ave(:,1),sub_timestamp),:)=[]; 
          elseif d(r)<0 %rew higher counts that ave
                sub_occ(r)=N_ave(r);%keep ave counts 
                %Subsampling rew spk 
                %randomly chose the same amount of timestamps in that bin
                %in rew than in ave
                rew_bin = find(pos_rew(:,3)==r);
                bin = pos_rew(rew_bin,:);
                k = sort(randperm(numel(rew_bin),N_rew(r)-N_ave(r)));
                sub_timestamp = bin(k,1); 
                %Remove those bin time stamp from  bspks_rew
                bspks_rew(ismember(bspks_rew(:,1),sub_timestamp),:)=[]; 
          elseif d(r)==0   % equal counts
                sub_occ(r)=N_ave(r);
                % do not touch spk  
          end
                   
      end 

                OccMap = sub_occ*dt;
                OccMap(OccMap==0) = eps;
                miSmooth = 9;OccMap = imgaussfilt(OccMap,2,'FilterSize',miSmooth,'Padding','circular');

                %Create a N_spike map
                tspk_tr(tspk_tr(:,2)==0,:)=[];
                Xs = []; Ys = [];
                for i = 1:size(tspk_tr,1)
                    x = pos_tr(pos_tr(:,1)== tspk_tr(i,1),2);
                    y = pos_tr(pos_tr(:,1)== tspk_tr(i,1),3);
    
                    if tspk_tr(i,2)== 1 
                        Xs = [Xs;x];
                        Ys = [Ys;y];
                    else
                        Xs = [Xs;repelem(x,tspk_tr(i,2))'];
                        Ys = [Ys;repelem(y,tspk_tr(i,2))'];
                    end    
                end 
                [Ntr] = histcounts2(Xs,Ys,Xedges,Yedges);
                Ntr = imgaussfilt(Ntr,2,'FilterSize',miSmooth,'Padding','replicate'); 
                Ntr(indx_mascara_tr) = nan; 
                
                tspk_ts(tspk_ts(:,2)==0,:)=[];
                Xs = []; Ys = [];
                for i = 1:size(tspk_ts,1)
                    x = pos_ts(pos_ts(:,1)== tspk_ts(i,1),2);
                    y = pos_ts(pos_ts(:,1)== tspk_ts(i,1),3);
    
                    if tspk_ts(i,2)== 1 
                        Xs = [Xs;x];
                        Ys = [Ys;y];
                    else
                        Xs = [Xs;repelem(x,tspk_ts(i,2))'];
                        Ys = [Ys;repelem(y,tspk_ts(i,2))'];
                    end    
                end 
                [Nts] = histcounts2(Xs,Ys,Xedges,Yedges);
                Nts = imgaussfilt(Nts,2,'FilterSize',miSmooth,'Padding','replicate');
                Nts(indx_mascara_ts) = nan;
                
                RateMapPROBE = Ntr./OccMap;
                RateMapTEST = Nts./OccMap;






end