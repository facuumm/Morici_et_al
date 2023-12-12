function [id , percentage , counts] = from_members_to_placecells(path)
% ----- Routine for percentage calculation of PC in assemblies members -----
% This function iterates inside the path you introduce and will upload the
% mat files containing dorsal and ventral hippocampal place-cells ids and
% will match those ones to the assemblies members ids defined.
%
%
% syntax: waveform_parameters(path,sf)
%
% INPUTS
% path: cell, contains the path of the sessions you want to analyze.
%       This scripts get inside each subfolder of the introduced path and
%       will get into the 'Spikesorting' folter where it will look for the
%       output of loadSpikes.m
%
% OUTPUS
% id: contains the place-cells ids organaized in the following
%           manner:   id.TYPE_OF_ASSEMBLIE.CONDITION.STRUCTURE_WHERE_PC_CAME
%
% percentage: contains the percentage of place-cells members relative to
%             the total members from the assembly. Same organization as it
%             described above.
%
% counts: contains the members number in each condition and structure
%
% Facundo Morici 12/2023


id.vHPC.aversive.vHPC = [];     percentage.vHPC.aversive.vHPC = [];
id.joint.aversive.vHPC = [];    percentage.joint.aversive.vHPC = [];
id.vHPC.reward.vHPC = [];       percentage.vHPC.reward.vHPC = [];
id.joint.reward.vHPC = [];      percentage.joint.reward.vHPC = [];
id.dHPC.aversive.dHPC = [];     percentage.dHPC.aversive.dHPC = [];
id.joint.aversive.dHPC = [];    percentage.joint.aversive.dHPC = [];
id.dHPC.reward.dHPC = [];       percentage.dHPC.reward.dHPC = [];
id.joint.reward.dHPC = [];      percentage.joint.reward.dHPC = [];
counts.aversive.dHPC = [];  counts.aversive.vHPC = [];  counts.aversive.joint.dHPC = [];
counts.reward.dHPC = [];    counts.reward.vHPC = [];    counts.reward.joint.dHPC = [];
counts.reward.joint.vHPC = [];  counts.aversive.joint.vHPC = [];

for tt = 1:length(path)
    %List of folders from the path
    files = dir(path{tt});
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
    clear files dirFlags
    
    for t = 1 : length(subFolders)-2
        session = [subFolders(t+2).folder,'\',subFolders(t+2).name,'\Spikesorting'];
        cd(session)
        
        %% load clusters used to detect assemblies
        x = dir([cd,'\clusters_included_in*.mat']);
        if not(isempty(x))
            load(x(end).name); clear x
        end
        %% load assemblies patterns
        x = dir([cd,'\dorsalventral_assemblies_reward*.mat']);
        if not(isempty(x))
            load(x(end).name)
            clear x
            TH.reward = Th;
            clear Th criteria_fr criteria_n pat x
        end
        clear x
        
        x = dir([cd,'\dorsalventral_assemblies_aversive*.mat']);
        if not(isempty(x))
            load(x(end).name)
            clear x
            TH.aversive = Th;
            clear Th criteria_fr criteria_n pat x
        end
        clear x
        if exist('TH')
            %% Definition of dHPC, vHPC or joint assemblies
            % Aversive
            cond1 =  sum(TH.aversive(1:size(clusters.dHPC,1),:),1)>0; %checking of dHPC SU
            cond2 =  sum(TH.aversive(size(clusters.dHPC,1)+1:end,:),1)>0; %checking of vHPC SU
            cond.dHPC.aversive = and(cond1 , not(cond2));
            cond.vHPC.aversive = and(cond2 , not(cond1));
            cond.both.aversive = and(cond1 , cond2); clear cond1 cond2
            
            TH.dHPC.aversive = TH.aversive(:,cond.dHPC.aversive);
            TH.vHPC.aversive = TH.aversive(:,cond.vHPC.aversive);
            TH.both.aversive = TH.aversive(:,cond.both.aversive);
            
            % Reward
            cond1 =  sum(TH.reward(1:size(clusters.dHPC,1),:),1)>0; %checking of dHPC SU
            cond2 =  sum(TH.reward(size(clusters.dHPC,1)+1:end,:),1)>0; %checking of vHPC SU
            cond.dHPC.reward = and(cond1 , not(cond2));
            cond.vHPC.reward = and(cond2 , not(cond1));
            cond.both.reward = and(cond1 , cond2); clear cond1 cond2
            
            TH.dHPC.reward = TH.reward(:,cond.dHPC.reward);
            TH.vHPC.reward = TH.reward(:,cond.vHPC.reward);
            TH.both.reward = TH.reward(:,cond.both.reward);
            
            clear cond
            
            %% load information of place-cell analysis
            if isfile('vHPC_pc.mat') %vHPC place-cells
                x = dir([cd,'\vHPC_pc.mat']);
                load(x(end).name)
                
                % vHPC assemblies Aversive
                for i = 1 : size(TH.vHPC.aversive,2)
                    count = 0;
                    c = clusters.all(TH.vHPC.aversive(:,i));
                    
                    for ii = 1:size(vHPC,2)
                        tmp = c == vHPC{ii}.id;
                        if sum(tmp)>0
                            count = count+1;
                            id.vHPC.aversive.vHPC = [id.vHPC.aversive.vHPC ; vHPC{ii}.id];
                        end
                    end
                    
                    percentage.vHPC.aversive.vHPC = [percentage.vHPC.aversive.vHPC ; (count/sum(TH.vHPC.aversive(:,i)))];
                    counts.aversive.vHPC = [counts.aversive.vHPC ; sum(TH.vHPC.aversive(:,i))];
                    clear count c
                end
                
                % Joint assemblies Aversive
                for i = 1 : size(TH.both.aversive,2)
                    count = 0;
                    c = clusters.all(TH.both.aversive(:,i));
                    
                    for ii = 1:size(vHPC,2)
                        tmp = c == vHPC{ii}.id;
                        if sum(tmp)>0
                            count = count+1;
                            id.joint.aversive.vHPC = [id.joint.aversive.vHPC ; vHPC{ii}.id];
                        end
                    end
                    percentage.joint.aversive.vHPC = [percentage.joint.aversive.vHPC ; (count/sum(TH.both.aversive(size(clusters.dHPC,1)+1:end,i)))];
                    counts.aversive.joint.vHPC = [counts.aversive.joint.vHPC ; sum(TH.both.aversive(size(clusters.dHPC,1)+1:end,i))];
                    clear count c
                end
                
                % vHPC assemblies Reward
                for i = 1 : size(TH.vHPC.reward,2)
                    count = 0;
                    c = clusters.all(TH.vHPC.reward(:,i));
                    
                    for ii = 1:size(vHPC,2)
                        tmp = c == vHPC{ii}.id;
                        if sum(tmp)>0
                            count = count+1;
                            id.vHPC.reward.vHPC = [id.vHPC.reward.vHPC ; vHPC{ii}.id];
                        end
                    end
                    percentage.vHPC.reward.vHPC = [percentage.vHPC.reward.vHPC ; (count/sum(TH.vHPC.reward(:,i)))];
                    counts.reward.vHPC = [counts.reward.vHPC ; sum(TH.vHPC.reward(:,i))];
                    clear count c
                end
                
                % Joint assemblies Reward
                for i = 1 : size(TH.both.reward,2)
                    count = 0;
                    c = clusters.all(TH.both.reward(:,i));
                    
                    for ii = 1:size(vHPC,2)
                        tmp = c == vHPC{ii}.id;
                        if sum(tmp)>0
                            count = count+1;
                            id.vHPC.reward.vHPC = [id.vHPC.reward.vHPC ; vHPC{ii}.id];
                        end
                    end
                    percentage.joint.reward.vHPC = [percentage.joint.reward.vHPC ; (count/sum(TH.both.reward(size(clusters.dHPC,1)+1:end,i)))];
                    counts.reward.joint.vHPC = [counts.reward.joint.vHPC ; sum(TH.both.reward(size(clusters.dHPC,1)+1:end,i))];
                    clear count c
                end
                
            end
            
            if isfile('dHPC_pc.mat') %dHPC place-cells
                x = dir([cd,'\dHPC_pc.mat']);
                load(x(end).name)
                
                % dHPC assemblies Aversive
                for i = 1 : size(TH.dHPC.aversive,2)
                    count = 0;
                    c = clusters.all(TH.dHPC.aversive(:,i));
                    
                    for ii = 1:size(dHPC,2)
                        tmp = c == dHPC{ii}.id;
                        if sum(tmp)>0
                            count = count+1;
                            id.dHPC.aversive.dHPC = [id.dHPC.aversive.dHPC ; dHPC{ii}.id];
                        end
                    end
                    percentage.dHPC.aversive.dHPC = [percentage.dHPC.aversive.dHPC ; (count/sum(TH.dHPC.aversive(:,i)))];
                    counts.aversive.dHPC = [counts.aversive.dHPC ; sum(TH.dHPC.aversive(:,i))];
                    clear count c
                end
                
                % Joint assemblies Aversive
                for i = 1 : size(TH.both.aversive,2)
                    count = 0;
                    c = clusters.all(TH.both.aversive(:,i));
                    
                    for ii = 1:size(dHPC,2)
                        tmp = c == dHPC{ii}.id;
                        if sum(tmp)>0
                            count = count+1;
                            id.joint.aversive.dHPC = [id.joint.aversive.dHPC ; dHPC{ii}.id];
                        end
                    end
                    percentage.joint.aversive.dHPC = [percentage.joint.aversive.dHPC ; (count/sum(TH.both.aversive(1:size(clusters.dHPC,1),i)))];
                    counts.aversive.joint.dHPC = [counts.aversive.joint.dHPC ; sum(TH.both.aversive(1:size(clusters.dHPC,1),i))];
                    clear count c
                end
                
                % dHPC assemblies Reward
                for i = 1 : size(TH.dHPC.reward,2)
                    count = 0;
                    c = clusters.all(TH.dHPC.reward(:,i));
                    
                    for ii = 1:size(dHPC,2)
                        tmp = c == dHPC{ii}.id;
                        if sum(tmp)>0
                            count = count+1;
                            id.dHPC.reward.dHPC = [id.dHPC.reward.dHPC ; dHPC{ii}.id];
                        end
                    end
                    percentage.dHPC.reward.dHPC = [percentage.dHPC.reward.dHPC ; (count/sum(TH.dHPC.reward(:,i)))];
                    counts.reward.dHPC = [counts.reward.dHPC ; sum(TH.dHPC.reward(:,i))];
                    clear count c
                end
                
                % Joint assemblies Reward
                for i = 1 : size(TH.both.reward,2)
                    count = 0;
                    c = clusters.all(TH.both.reward(:,i));
                    
                    for ii = 1:size(dHPC,2)
                        tmp = c == dHPC{ii}.id;
                        if sum(tmp)>0
                            count = count+1;
                            id.joint.reward.vHPC = [id.joint.reward.dHPC ; dHPC{ii}.id];
                        end
                    end
                    percentage.joint.reward.dHPC = [percentage.joint.reward.dHPC ; (count/sum(TH.both.reward(1:size(clusters.dHPC,1),i)))];
                    counts.reward.joint.dHPC = [counts.reward.joint.dHPC ; sum(TH.both.reward(1:size(clusters.dHPC,1),i))];
                    clear count c
                end
            end
            
            clear dHPC vHPC i ii x clusters TH tmp
        end
    end

end
end
