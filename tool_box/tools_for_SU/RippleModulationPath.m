function RippleModulationPath(path)
% This function determine if a SU is modulated by the ripples.
%
% syntax = RippleModulationPath(path)
%
% INPUTS
% path: cell, contains the paths of the sessions you want to analyze.
%
% OUTPUT
% pInc = probability to be up-modulated
% pDec = probability to be down-modulated
% suprise = positive if pInc > pDec
%
% It save in each Spikesorting folder a mat file called 'RippleModulatedSU'
% It contains three structures with the following organization.
%
% pInc.dHPC
%     .vHPC
%     .dvHPC.dHPC
%           .vHPC
% pDec ...
% surp ...
%
% The first column contains the cluster ID
% The second define the cellulartype (1: Pyr, 0: Int)
% The third column, in pInc or pDec if is 1, then is up or down-modulated.
%                 , in surp it contains surprise value  
%
% other functions: poissonTest (Eran Stark)
%                  Restrict and SubsSubtractIntervals (FMAtoolbox)
%
% Morci Juan Facundo 01/2024





%% Main loop, to iterate across sessions
for tt = 1:length(path)
    %List of folders from the path
    files = dir(path{tt});
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
    clear files dirFlags
    num_assembliesR = [];
    num_assembliesA = [];
    for t = 1 : length(subFolders)-2
        disp(['-- Initiating analysis of folder #' , num2str(t) , ' from rat #',num2str(tt) , ' --'])
        session = [subFolders(t+2).folder,'\',subFolders(t+2).name];
        cd(session)
        %Loading TS of the sessions
        disp('Uploading session time stamps')
        x = dir([cd,'\*.cat.evt']);
        segments = readtable([cd,'\',x.name],'FileType','text');
        clear x
        
        %% load sleep states
        disp('Uploading sleep scoring')
        x = dir([cd,'\*-states.mat']);    states = load([cd,'\',x.name]);    states = states.states;
        REM.all = ToIntervals(states==5);    NREM.all = ToIntervals(states==3);    WAKE.all = ToIntervals(states==1);
        %         REM.all(REM.all(:,2)-REM.all(:,1)<60,:) = [];
        clear x states
        
        %% Load ripples
        if isfile('ripplesD_customized2.csv')
            ripplesD = table2array(readtable('ripplesD_customized2.csv'));
            condD = true;
        else
            condD = false;
        end
        
        if isfile('ripplesV_customized2.csv')
            ripplesV = table2array(readtable('ripplesV_customized2.csv'));
            condV = true;
        else
            condV = false;
        end
        
        
        %% Spikes
        % Load Units
        disp('Uploading Spiking activity')
        cd 'Spikesorting'
        Spks = double([readNPY('spike_clusters.npy') readNPY('spike_times.npy')]);
        K = tdfread('cluster_group.tsv'); % import clusters ID and groups (MUA,Noise,Good)
        Kinfo = tdfread('cluster_info.tsv'); % import information of clusters
        K = [K.cluster_id(K.group(:,1) == 'g') , Kinfo.ch(K.group(:,1) == 'g'),]; % defining good clusters
        % Load neuronal classification
        load('Cell_type_classification')
        K = [K , Cell_type_classification(:,6:7)];
        group_dHPC = K(K(:,2) > 63,:);
        group_vHPC = K(K(:,2) <= 63,:);
        
        %Loop to select dorsal or ventral LFP and SU
        % z=1 --> dorsal
        % z=2 --> ventral
        for z = 1:2
            if z == 1
                spks_dHPC = Spks(ismember(Spks(:,1),group_dHPC(:,1)),:); %keep spks from good clusters
            else
                spks_vHPC = Spks(ismember(Spks(:,1),group_vHPC(:,1)),:); %keep spks from good clusters
            end
        end
        clear z
        spks_vHPC(:,2) = double(spks_vHPC(:,2))./20000;
        spks_dHPC(:,2) = double(spks_dHPC(:,2))./20000;
        Spks(:,2) = double(Spks(:,2))./20000;
        
        
        cellulartype = [K(:,1) , K(:,4)]; % 1 if is Pyr and 0 if is Int
        
        
        %% Counting the Number f SU
        clusters.dHPC = [];
        for ii=1:size(group_dHPC,1)
            cluster = group_dHPC(ii,1);
            Cellulartype = logical(cellulartype(cellulartype(:,1) == cluster,2));
            if Cellulartype
                clusters.dHPC = [clusters.dHPC ; cluster , true];
            else
                clusters.dHPC = [clusters.dHPC ; cluster , false];
            end
            clear tmp cluster Cellulartype fr1 fr2 fr3 fr4 fr5 r a
        end
        
        clusters.vHPC = [];
        for ii=1:size(group_vHPC,1)
            cluster = group_vHPC(ii,1);
            Cellulartype = logical(cellulartype(cellulartype(:,1) == cluster,2));
            if Cellulartype
                clusters.vHPC = [clusters.vHPC ; cluster , true];
            else
                clusters.vHPC = [clusters.vHPC ; cluster , false];
            end
            clear tmp cluster Cellulartype fr1 fr2 fr3 fr4 fr5 r a
        end
        
        
        %% PoisonTest
        
        if condD
            pInc.dHPC = [];
            pDec.dHPC = [];
            surp.dHPC = [];
            iii = Restrict(ripplesD,NREM.all);
            bufferedripples=[iii(:,1)-0.1 iii(:,3)+0.1];
            [baseline,ind]=SubtractIntervals(NREM.all,bufferedripples,'strict','on');
            for i = 1: size(clusters.dHPC,1)
                cluster = clusters.dHPC(i,:);
                spks=Restrict(Spks(Spks(:,1)==cluster(1),2),NREM.all);
                
                totalbaselinetime=sum(baseline(:,2)-baseline(:,1));
                baselinespikes=Restrict(spks(:,1),baseline);
                
                % Restrict further to the third of ripples with the highest amplitude
                totalrippletime=sum(iii(:,3)-iii(:,1));
                ripplespikes=Restrict(spks(:,1),[iii(:,1) iii(:,3)]);
                ncellbaselinespikes=length(baselinespikes);
                ncellripplespikes=length(ripplespikes);
                
                if ncellbaselinespikes~=0 & ncellripplespikes~=0
                    [Inc Dec sur] = poissonTest(ncellbaselinespikes/totalbaselinetime,ncellripplespikes,totalrippletime);
                    pInc.dHPC = [pInc.dHPC ; cluster , Inc<0.01/2];
                    pDec.dHPC = [pDec.dHPC ; cluster , Dec<0.01/2];
                    surp.dHPC = [surp.dHPC ; cluster , sur];
                    clear cluster Inc Dec sur
                else
                    [Inc Dec sur] = poissonTest(ncellbaselinespikes/totalbaselinetime,ncellripplespikes,totalrippletime);
                    pInc.dHPC = [pInc.dHPC ; cluster , NaN];
                    pDec.dHPC = [pDec.dHPC ; cluster , NaN];
                    surp.dHPC = [surp.dHPC ; cluster , NaN];
                    clear cluster Inc Dec sur
                end
                clear totalrippletime ripplespikes ncellbaselinespikes ncellripplespikes
                clear baselinespikes totalbaselinetime
            end
            clear bufferedripples iii
        else
            pInc.dHPC = [];
            pDec.dHPC = [];
            surp.dHPC = [];
            clear totalrippletime ripplespikes ncellbaselinespikes ncellripplespikes
            clear baselinespikes totalbaselinetime bufferedripples iii
        end
        
        
        if condV
            pInc.vHPC = [];
            pDec.vHPC = [];
            surp.vHPC = [];
            iii = Restrict(ripplesV,NREM.all);
            bufferedripples=[iii(:,1)-0.1 iii(:,3)+0.1];
            [baseline,ind]=SubtractIntervals(NREM.all,bufferedripples,'strict','on');
            for i = 1: size(clusters.vHPC,1)
                cluster = clusters.vHPC(i,:);
                spks=Restrict(Spks(Spks(:,1)==cluster(1),2),NREM.all);
                
                totalbaselinetime=sum(baseline(:,2)-baseline(:,1));
                baselinespikes=Restrict(spks(:,1),baseline);
                
                % Restrict further to the third of ripples with the highest amplitude
                totalrippletime=sum(iii(:,3)-iii(:,1));
                ripplespikes=Restrict(spks(:,1),[iii(:,1) iii(:,3)]);
                ncellbaselinespikes=length(baselinespikes);
                ncellripplespikes=length(ripplespikes);
                
                if ncellbaselinespikes~=0 & ncellripplespikes~=0
                    [Inc Dec sur] = poissonTest(ncellbaselinespikes/totalbaselinetime,ncellripplespikes,totalrippletime);
                    pInc.vHPC = [pInc.vHPC ; cluster , Inc<0.01/2];
                    pDec.vHPC = [pDec.vHPC ; cluster , Dec<0.01/2];
                    surp.vHPC = [surp.vHPC ; cluster , sur];
                    clear cluster Inc Dec sur
                else
                    [Inc Dec sur] = poissonTest(ncellbaselinespikes/totalbaselinetime,ncellripplespikes,totalrippletime);
                    pInc.vHPC = [pInc.vHPC ; cluster , NaN];
                    pDec.vHPC = [pDec.vHPC ; cluster , NaN];
                    surp.vHPC = [surp.vHPC ; cluster , NaN];
                    clear cluster Inc Dec sur
                end
                clear totalrippletime ripplespikes ncellbaselinespikes ncellripplespikes
                clear baselinespikes totalbaselinetime
            end
            clear bufferedripples iii
        else
            pInc.vHPC = [];
            pDec.vHPC = [];
            surp.vHPC = [];
            clear totalrippletime ripplespikes ncellbaselinespikes ncellripplespikes
            clear baselinespikes totalbaselinetime bufferedripples iii
        end
        
        
        if and(condD,condV)
            pInc.dvHPC.dHPC = [];
            pDec.dvHPC.dHPC = [];
            surp.dvHPC.dHPC = [];
            iii = Restrict(ripplesV,NREM.all);
            bufferedripples=[iii(:,1)-0.1 iii(:,3)+0.1];
            [baseline,ind]=SubtractIntervals(NREM.all,bufferedripples,'strict','on');
            for i = 1: size(clusters.dHPC,1)
                cluster = clusters.dHPC(i,:);
                spks=Restrict(Spks(Spks(:,1)==cluster(1),2),NREM.all);
                
                totalbaselinetime=sum(baseline(:,2)-baseline(:,1));
                baselinespikes=Restrict(spks(:,1),baseline);
                
                % Restrict further to the third of ripples with the highest amplitude
                totalrippletime=sum(iii(:,3)-iii(:,1));
                ripplespikes=Restrict(spks(:,1),[iii(:,1) iii(:,3)]);
                ncellbaselinespikes=length(baselinespikes);
                ncellripplespikes=length(ripplespikes);
                
                if ncellbaselinespikes~=0 & ncellripplespikes~=0
                    [Inc Dec sur] = poissonTest(ncellbaselinespikes/totalbaselinetime,ncellripplespikes,totalrippletime);
                    pInc.dvHPC.dHPC = [pInc.dvHPC.dHPC ; cluster , Inc<0.01/2];
                    pDec.dvHPC.dHPC = [pDec.dvHPC.dHPC ; cluster , Dec<0.01/2];
                    surp.dvHPC.dHPC = [surp.dvHPC.dHPC ; cluster , sur];
                    clear cluster Inc Dec sur
                else
                    [Inc Dec sur] = poissonTest(ncellbaselinespikes/totalbaselinetime,ncellripplespikes,totalrippletime);
                    pInc.dvHPC.dHPC = [pInc.dvHPC.dHPC ; cluster , NaN];
                    pDec.dvHPC.dHPC = [pDec.dvHPC.dHPC ; cluster , NaN];
                    surp.dvHPC.dHPC = [surp.dvHPC.dHPC ; cluster , NaN];
                    clear cluster Inc Dec sur
                end
                clear totalrippletime ripplespikes ncellbaselinespikes ncellripplespikes
                clear baselinespikes totalbaselinetime
            end
            clear bufferedripples iii
            
            % vHPC SU
            pInc.dvHPC.vHPC = [];
            pDec.dvHPC.vHPC = [];
            surp.dvHPC.vHPC = [];
            iii = Restrict(ripplesD,NREM.all);
            bufferedripples=[iii(:,1)-0.1 iii(:,3)+0.1];
            [baseline,ind]=SubtractIntervals(NREM.all,bufferedripples,'strict','on');
            for i = 1: size(clusters.vHPC,1)
                cluster = clusters.vHPC(i,:);
                spks=Restrict(Spks(Spks(:,1)==cluster(1),2),NREM.all);
                
                totalbaselinetime=sum(baseline(:,2)-baseline(:,1));
                baselinespikes=Restrict(spks(:,1),baseline);
                
                % Restrict further to the third of ripples with the highest amplitude
                totalrippletime=sum(iii(:,3)-iii(:,1));
                ripplespikes=Restrict(spks(:,1),[iii(:,1) iii(:,3)]);
                ncellbaselinespikes=length(baselinespikes);
                ncellripplespikes=length(ripplespikes);
                
                if ncellbaselinespikes~=0 & ncellripplespikes~=0
                    [Inc Dec sur] = poissonTest(ncellbaselinespikes/totalbaselinetime,ncellripplespikes,totalrippletime);
                    pInc.dvHPC.vHPC = [pInc.dvHPC.vHPC ; cluster , Inc<0.01/2];
                    pDec.dvHPC.vHPC = [pDec.dvHPC.vHPC ; cluster , Dec<0.01/2];
                    surp.dvHPC.vHPC = [surp.dvHPC.vHPC ; cluster , sur];
                    clear cluster Inc Dec sur
                else
                    [Inc Dec sur] = poissonTest(ncellbaselinespikes/totalbaselinetime,ncellripplespikes,totalrippletime);
                    pInc.dvHPC.vHPC = [pInc.dvHPC.vHPC ; cluster , NaN];
                    pDec.dvHPC.vHPC = [pDec.dvHPC.vHPC ; cluster , NaN];
                    surp.dvHPC.vHPC = [surp.dvHPC.vHPC ; cluster , NaN];
                    clear cluster Inc Dec sur
                end
                clear totalrippletime ripplespikes ncellbaselinespikes ncellripplespikes
                clear baselinespikes totalbaselinetime
            end
            clear bufferedripples iii
            
            
        else
            pInc.dvHPC.dHPC = [];
            pDec.dvHPC.dHPC = [];
            surp.dvHPC.dHPC = [];
            
            pInc.dvHPC.vHPC = [];
            pDec.dvHPC.vHPC = [];
            surp.dvHPC.vHPC = [];
            
            clear totalrippletime ripplespikes ncellbaselinespikes ncellripplespikes
            clear baselinespikes totalbaselinetime bufferedripples iii
        end
        
        save([cd,'\RippleModulatedSU.mat'], 'pInc', 'pDec' , 'surp');
        clear pInc pDec surp Spks ripplesD ripplesV cellulartype Cellulartype
        clear condV condD group_dHPC group_vHPC K Kinfo NREM WAKE REM
        clear spks_dHPC spks_vHPC ii i Cell_type_classification clusters segments spks
    end
end
end