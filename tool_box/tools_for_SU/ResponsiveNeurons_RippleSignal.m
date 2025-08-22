function [dorsal , ventral] = ResponsiveNeurons_RippleSignal(path)
% This function iterates inside the path uploading the Neurons, their
% responsivnes to shocks and the SWR polarity (in SWR_sign.csv).
%
% INPUTS
% path: cell, inside each cell you should put the path of each session you
%       want to analyze
%
% OUTPUT
% dorsal, ventral: matrix, contains the id, shock_responsivnes, channel
%                  where they was detected and the SWR polarity of that
%                  channel.
%
% Morci Juan Facundo 08/2025

% Variables to use in the script
criteria_fr = 0;
criteria_type = 0; % criteria for celltype (0:pyr, 1:int, 2:all)

% Storage variables
dorsal.shock = [];    ventral.shock = [];
dorsal.valve = [];    ventral.valve = [];
%% Main loop, to iterate across sessions
for tt = 1:length(path)
    %List of folders from the path
    files = dir(path{tt});
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
    clear files dirFlags
    
    for t = 1 : length(subFolders)-2
        disp(['-- Initiating analysis of folder #' , num2str(t) , ' from rat #',num2str(tt) , ' --'])
        session = [subFolders(t+2).folder,'\',subFolders(t+2).name];
        cd(session)
        
        %% Load Spiking Groups
        disp('Loading Spiking Groups')
        x = dir([cd,'\',subFolders(t+2).name,'.xml']);
        SpkGrps = LoadXml([x.folder,'\',x.name]); clear x
        SpkGrps = SpkGrps.SpkGrps;
        
        %% Load Session organization
        disp('Loading TimeStamps Session Structure')
        load('session_organization.mat')
        
        %% Load Ripple Polarity
        cd 'Spikesorting'
        disp('Loading Ripples Polarity')
        RipplePolarity = readmatrix([cd,'\SWR_sign.csv']);
        
        %% Loading Channel Positions
        disp('Loading Channel Positions')
        x1 = dir([cd,'\channel_map.npy']);
        x2 = dir([cd,'\channel_positions.npy']);
        
        x1 = readNPY([x1.folder,'\',x1.name]);
        x2 = readNPY([x2.folder,'\',x2.name]);
        
        ChannelPositions = [x1' , x2]; clear x1 x2
        
        spike_clusters = readNPY([cd,'\spike_clusters.npy']);
        spike_templates = readNPY([cd,'\spike_templates.npy']);
        spike_templates = spike_templates + 1;
        templates = readNPY([cd,'\templates.npy']);
        
        
        %% Spikes
        % Load Units
        disp('Uploading Spiking activity')        
        [clusters , numberD , numberV , spks , spks_dHPC , spks_vHPC , cellulartype] = load_SU_FM(cd,criteria_type,criteria_fr,aversiveTS_run./1000,rewardTS_run./1000);
        
        %% Shock
        if exist('dHPC_Shock_VF.mat')
            load('dHPC_Shock_VF.mat')
            for i = 1 : size(dHPC_resp.id,1)
                cluster = dHPC_resp.id(i);
                tmp = cellulartype(:,1) == cluster;
                ch = cellulartype(tmp,3);
                Polarity = RipplePolarity(RipplePolarity(:,1) == ch , 2);
                depth = ChannelPositions(ChannelPositions(:,1) == ch,3);
                
                % Detection of Spiking group
                for ii = 1:length(SpkGrps)
                    if any(SpkGrps(ii).Channels == ch)
                        Channels = SpkGrps(ii).Channels;
                        D = [];
                        for iii = 1 : length(Channels)
                            D = [D ; ChannelPositions(ChannelPositions(:,1) == Channels(iii),2:3)];
                        end
                        groupChannels = [Channels'+1 , D]; clear Channels D
                        break
                    end
                end
                
                % Calculation of triangulization
                idx = spike_clusters == cluster;
                templatesClu = unique(spike_templates(idx)); clear idx
                
                FinalTemplate = templates(templatesClu , : , :); clear templatesClu
                FinalTemplate = FinalTemplate(:,:,groupChannels(:,1));
                
                [depth_um, depth_norm] = getUnitDepth(FinalTemplate, groupChannels);
                
                dorsal.shock = [dorsal.shock ; double(dHPC_resp.id(i)) , double(dHPC_resp.resp_ave(i)) , double(ch) , double(Polarity) , double(depth) , double(depth_norm)];
                clear tmp ch Polarity depth FinalTemplate groupChannels depth_um cluster
            end
        end
        
        if exist('vHPC_Shock_VF.mat')
            load('vHPC_Shock_VF.mat')
            for i = 1 : size(vHPC_resp.id,1)
                cluster = vHPC_resp.id(i);
                tmp = cellulartype(:,1) == cluster;
                ch = cellulartype(tmp,3);
                Polarity = RipplePolarity(RipplePolarity(:,1) == ch , 2);
                depth = ChannelPositions(ChannelPositions(:,1) == ch,3);
                
                % Detection of Spiking group
                for ii = 1:length(SpkGrps)
                    if any(SpkGrps(ii).Channels == ch)
                        Channels = SpkGrps(ii).Channels;
                        D = [];
                        for iii = 1 : length(Channels)
                            D = [D ; ChannelPositions(ChannelPositions(:,1) == Channels(iii),2:3)];
                        end
                        groupChannels = [Channels'+1 , D]; clear Channels D
                        break
                    end
                end                
                % Calculation of triangulization
                idx = spike_clusters == cluster;
                templatesClu = unique(spike_templates(idx)); clear idx
                
                FinalTemplate = templates(templatesClu , : , :); clear templatesClu
                FinalTemplate = FinalTemplate(:,:,groupChannels(:,1));
                
                [depth_um, depth_norm] = getUnitDepth(FinalTemplate, groupChannels);                
                
                ventral.shock = [ventral.shock ; double(vHPC_resp.id(i)) , double(vHPC_resp.resp_ave(i)) , double(ch) , double(Polarity) , double(depth) , double(depth_norm)];
                clear tmp ch Polarity depth depth_um depth_norm cluster
            end
        end
        
        
        
        %% Valve
        if exist('dHPC_valve.mat')
            load('dHPC_valve.mat')
            for i = 1 : size(dHPC_valve.id,1)
                cluster = dHPC_valve.id(i);
                tmp = cellulartype(:,1) == cluster;
                ch = cellulartype(tmp,3);
                Polarity = RipplePolarity(RipplePolarity(:,1) == ch , 2);
                depth = ChannelPositions(ChannelPositions(:,1) == ch,3);
                
                % Detection of Spiking group
                for ii = 1:length(SpkGrps)
                    if any(SpkGrps(ii).Channels == ch)
                        Channels = SpkGrps(ii).Channels;
                        D = [];
                        for iii = 1 : length(Channels)
                            D = [D ; ChannelPositions(ChannelPositions(:,1) == Channels(iii),2:3)];
                        end
                        groupChannels = [Channels'+1 , D]; clear Channels D
                        break
                    end
                end
                
                % Calculation of triangulization
                idx = spike_clusters == cluster;
                templatesClu = unique(spike_templates(idx)); clear idx
                
                FinalTemplate = templates(templatesClu , : , :); clear templatesClu
                FinalTemplate = FinalTemplate(:,:,groupChannels(:,1));
                
                [depth_um, depth_norm] = getUnitDepth(FinalTemplate, groupChannels);
                
                dorsal.valve = [dorsal.valve ; double(dHPC_valve.id(i)) , double(dHPC_valve.responssiveness(i)) , double(ch) , double(Polarity) , double(depth) , double(depth_norm)];
                clear tmp ch Polarity depth FinalTemplate groupChannels depth_um cluster
            end
        end
        
        if exist('vHPC_valve.mat')
            load('vHPC_valve.mat')
            for i = 1 : size(vHPC_valve.id,1)
                cluster = vHPC_valve.id(i);
                tmp = cellulartype(:,1) == cluster;
                ch = cellulartype(tmp,3);
                Polarity = RipplePolarity(RipplePolarity(:,1) == ch , 2);
                depth = ChannelPositions(ChannelPositions(:,1) == ch,3);
                
                % Detection of Spiking group
                for ii = 1:length(SpkGrps)
                    if any(SpkGrps(ii).Channels == ch)
                        Channels = SpkGrps(ii).Channels;
                        D = [];
                        for iii = 1 : length(Channels)
                            D = [D ; ChannelPositions(ChannelPositions(:,1) == Channels(iii),2:3)];
                        end
                        groupChannels = [Channels'+1 , D]; clear Channels D
                        break
                    end
                end                
                % Calculation of triangulization
                idx = spike_clusters == cluster;
                templatesClu = unique(spike_templates(idx)); clear idx
                
                FinalTemplate = templates(templatesClu , : , :); clear templatesClu
                FinalTemplate = FinalTemplate(:,:,groupChannels(:,1));
                
                [depth_um, depth_norm] = getUnitDepth(FinalTemplate, groupChannels);                
                
                ventral.valve = [ventral.valve ; double(vHPC_valve.id(i)) , double(vHPC_valve.responssiveness(i)) , double(ch) , double(Polarity) , double(depth) , double(depth_norm)];
                clear tmp ch Polarity depth depth_um depth_norm cluster
            end
        end        
        
        
        
    end
    disp(' ')
    clear A aversiveTS aversiveTS_run baselineTS rewardTS rewardTS_run
    clear is K Kinfo group_dHPC group_vHPC matrixC matrixD matrixV
    clear spiketrains_dHPC spiketrains_vHPC opts MUA
    clear patterns Thresholded i  ii numberD numberV movement cross crossN
    clear Spikes bins Clusters Shocks_filt Rewards_filt config n_SU_D n_SU_V
    clear clusters coordinated coordinated_ripple_bursts coordinatedV
    clear cooridnated_event coordinatedV_refined coordinatedV_refined
    clear ripple_bursts ripple_event ripplesD ripplesV
    clear spks spks_dHPC spks_vHPC ripples cooridnated_event
    clear cooridnated_eventDV cooridnated_eventVD segments movement RD RV RB
    clear dHPC vHPC1 vHPC2 vHPC3 vHPC4 vHPC5 distribution
end

%% ==========================
%  Shock
%  ==========================

% dHPC counts
counts.dHPC.up.shock   = sum(dorsal.shock(:,4) == 1 & dorsal.shock(:,2) == 1);
counts.dHPC.up.no      = sum(dorsal.shock(:,4) == 1 & dorsal.shock(:,2) ~= 1);
counts.dHPC.down.shock = sum(dorsal.shock(:,4) == 2 & dorsal.shock(:,2) == 1);
counts.dHPC.down.no    = sum(dorsal.shock(:,4) == 2 & dorsal.shock(:,2) ~= 1);

totalUp_d   = sum(dorsal.shock(:,4) == 1);
totalDown_d = sum(dorsal.shock(:,4) == 2);

percentage.dHPC.up.shock   = counts.dHPC.up.shock   / totalUp_d   * 100;
percentage.dHPC.up.no      = counts.dHPC.up.no      / totalUp_d   * 100;
percentage.dHPC.down.shock = counts.dHPC.down.shock / totalDown_d * 100;
percentage.dHPC.down.no    = counts.dHPC.down.no    / totalDown_d * 100;

% vHPC counts
counts.vHPC.up.shock   = sum(ventral.shock(:,4) == 1 & ventral.shock(:,2) == 1);
counts.vHPC.up.no      = sum(ventral.shock(:,4) == 1 & ventral.shock(:,2) ~= 1);
counts.vHPC.down.shock = sum(ventral.shock(:,4) == 2 & ventral.shock(:,2) == 1);
counts.vHPC.down.no    = sum(ventral.shock(:,4) == 2 & ventral.shock(:,2) ~= 1);

totalUp_v   = sum(ventral.shock(:,4) == 1);
totalDown_v = sum(ventral.shock(:,4) == 2);

percentage.vHPC.up.shock   = counts.vHPC.up.shock   / totalUp_v   * 100;
percentage.vHPC.up.no      = counts.vHPC.up.no      / totalUp_v   * 100;
percentage.vHPC.down.shock = counts.vHPC.down.shock / totalDown_v * 100;
percentage.vHPC.down.no    = counts.vHPC.down.no    / totalDown_v * 100;

% Plot dHPC
figure('Name','dHPC');

% --- Up
subplot(1,2,1);
vals = [percentage.dHPC.up.shock, percentage.dHPC.up.no];
nums = [counts.dHPC.up.shock, counts.dHPC.up.no];
b = bar(vals);
title('Up');
ylim([0 100]);
set(gca,'XTickLabel',{'Shock','No-Shock'});
ylabel('%');
for i = 1:numel(vals)
    text(i, vals(i)+3, sprintf('%.1f%% (%d/%d)', vals(i), nums(i), sum(dorsal(:,4)==1)), ...
        'HorizontalAlignment','center','FontSize',8);
end

% --- Down
subplot(1,2,2);
vals = [percentage.dHPC.down.shock, percentage.dHPC.down.no];
nums = [counts.dHPC.down.shock, counts.dHPC.down.no];
b = bar(vals);
title('Down');
ylim([0 100]);
set(gca,'XTickLabel',{'Shock','No-Shock'});
ylabel('%');
for i = 1:numel(vals)
    text(i, vals(i)+3, sprintf('%.1f%% (%d/%d)', vals(i), nums(i), sum(dorsal(:,4)==2)), ...
        'HorizontalAlignment','center','FontSize',8);
end

% Plot vHPC
figure('Name','vHPC');

% --- Up
subplot(1,2,1);
vals = [percentage.vHPC.up.shock, percentage.vHPC.up.no];
nums = [counts.vHPC.up.shock, counts.vHPC.up.no];
b = bar(vals);
title('Up');
ylim([0 100]);
set(gca,'XTickLabel',{'Shock','No-Shock'});
ylabel('%');
for i = 1:numel(vals)
    text(i, vals(i)+3, sprintf('%.1f%% (%d/%d)', vals(i), nums(i), sum(ventral(:,4)==1)), ...
        'HorizontalAlignment','center','FontSize',8);
end

% --- Down
subplot(1,2,2);
vals = [percentage.vHPC.down.shock, percentage.vHPC.down.no];
nums = [counts.vHPC.down.shock, counts.vHPC.down.no];
b = bar(vals);
title('Down');
ylim([0 100]);
set(gca,'XTickLabel',{'Shock','No-Shock'});
ylabel('%');
for i = 1:numel(vals)
    text(i, vals(i)+3, sprintf('%.1f%% (%d/%d)', vals(i), nums(i), sum(ventral(:,4)==2)), ...
        'HorizontalAlignment','center','FontSize',8);
end

% Calculation of deepness
% dHPC ratio: -20 (deep) --> 1, -160 (superficial) --> 0
sitesY_Shock = double.shock(dorsal.shock(:,6));

ratioD.all = sitesY_Shock;
ratioD.Shock = ratioD.all(dorsal.shock(:,2)==1);
ratioD.noShock = ratioD.all(~(dorsal.shock(:,2)==1));

% vHPC ratio: -160 (deep) --> 1, -20 (superficial) --> 0
sitesY_Shock = 1-double(ventral.shock(:,6));% Invert so that deep = 1

ratioV.all = sitesY_Shock;  % deep = 1
ratioV.all = 1 - ratioV.all;
ratioV.Shock = ratioV.all(ventral.shock(:,2)==1);
ratioV.noShock = ratioV.all(~(ventral.shock(:,2)==1));


x = [ones(size(ratioV.noShock)) ; ones(size(ratioV.Shock))*2];
y = [ratioV.noShock ; ratioV.Shock];
figure,
cdfplot(ratioV.noShock),hold on
cdfplot(ratioV.Shock),hold on
% boxplot(y,x)
[h p] = ranksum(ratioV.noShock , ratioV.Shock)


x = [ones(size(ratioD.noShock)) ; ones(size(ratioD.Shock))*2];
y = [ratioD.noShock ; ratioD.Shock];
figure,
cdfplot(ratioD.noShock),hold on
cdfplot(ratioD.Shock),hold on
% boxplot(y,x)
[h p] = ranksum(ratioD.noShock , ratioD.Shock)



%% ==========================
%  Valve
%  ==========================

% dHPC counts
counts.dHPC.up.valve   = sum(dorsal.valve(:,4) == 1 & dorsal.valve(:,2) == 1);
counts.dHPC.up.no      = sum(dorsal.valve(:,4) == 1 & dorsal.valve(:,2) ~= 1);
counts.dHPC.down.valve = sum(dorsal.valve(:,4) == 2 & dorsal.valve(:,2) == 1);
counts.dHPC.down.no    = sum(dorsal.valve(:,4) == 2 & dorsal.valve(:,2) ~= 1);

totalUp_d   = sum(dorsal.shock(:,4) == 1);
totalDown_d = sum(dorsal.shock(:,4) == 2);

percentage.dHPC.up.valve   = counts.dHPC.up.valve   / totalUp_d   * 100;
percentage.dHPC.up.no      = counts.dHPC.up.no      / totalUp_d   * 100;
percentage.dHPC.down.valve = counts.dHPC.down.valve / totalDown_d * 100;
percentage.dHPC.down.no    = counts.dHPC.down.no    / totalDown_d * 100;

% vHPC counts
counts.vHPC.up.valve   = sum(ventral.valve(:,4) == 1 & ventral.valve(:,2) == 1);
counts.vHPC.up.no      = sum(ventral.valve(:,4) == 1 & ventral.valve(:,2) ~= 1);
counts.vHPC.down.valve = sum(ventral.valve(:,4) == 2 & ventral.valve(:,2) == 1);
counts.vHPC.down.no    = sum(ventral.valve(:,4) == 2 & ventral.valve(:,2) ~= 1);

totalUp_v   = sum(ventral.valve(:,4) == 1);
totalDown_v = sum(ventral.valve(:,4) == 2);

percentage.vHPC.up.valve   = counts.vHPC.up.valve   / totalUp_v   * 100;
percentage.vHPC.up.no      = counts.vHPC.up.no      / totalUp_v   * 100;
percentage.vHPC.down.valve = counts.vHPC.down.valve / totalDown_v * 100;
percentage.vHPC.down.no    = counts.vHPC.down.no    / totalDown_v * 100;

% Plot dHPC
figure('Name','dHPC');

% --- Up
subplot(1,2,1);
vals = [percentage.dHPC.up.valve, percentage.dHPC.up.no];
nums = [counts.dHPC.up.valve, counts.dHPC.up.no];
b = bar(vals);
title('Up');
ylim([0 100]);
set(gca,'XTickLabel',{'valve','No-valve'});
ylabel('%');
for i = 1:numel(vals)
    text(i, vals(i)+3, sprintf('%.1f%% (%d/%d)', vals(i), nums(i), sum(dorsal.valve(:,4)==1)), ...
        'HorizontalAlignment','center','FontSize',8);
end

% --- Down
subplot(1,2,2);
vals = [percentage.dHPC.down.valve, percentage.dHPC.down.no];
nums = [counts.dHPC.down.valve, counts.dHPC.down.no];
b = bar(vals);
title('Down');
ylim([0 100]);
set(gca,'XTickLabel',{'valve','No-valve'});
ylabel('%');
for i = 1:numel(vals)
    text(i, vals(i)+3, sprintf('%.1f%% (%d/%d)', vals(i), nums(i), sum(dorsal.valve(:,4)==2)), ...
        'HorizontalAlignment','center','FontSize',8);
end

% Plot vHPC
figure('Name','vHPC');

% --- Up
subplot(1,2,1);
vals = [percentage.vHPC.up.valve, percentage.vHPC.up.no];
nums = [counts.vHPC.up.valve, counts.vHPC.up.no];
b = bar(vals);
title('Up');
ylim([0 100]);
set(gca,'XTickLabel',{'valve','No-valve'});
ylabel('%');
for i = 1:numel(vals)
    text(i, vals(i)+3, sprintf('%.1f%% (%d/%d)', vals(i), nums(i), sum(ventral.valve(:,4)==1)), ...
        'HorizontalAlignment','center','FontSize',8);
end

% --- Down
subplot(1,2,2);
vals = [percentage.vHPC.down.valve, percentage.vHPC.down.no];
nums = [counts.vHPC.down.valve, counts.vHPC.down.no];
b = bar(vals);
title('Down');
ylim([0 100]);
set(gca,'XTickLabel',{'valve','No-valve'});
ylabel('%');
for i = 1:numel(vals)
    text(i, vals(i)+3, sprintf('%.1f%% (%d/%d)', vals(i), nums(i), sum(ventral.valve(:,4)==2)), ...
        'HorizontalAlignment','center','FontSize',8);
end

% Calculation of deepness
% dHPC ratio: -20 (deep) --> 1, -160 (superficial) --> 0
sitesY_valve = double(dorsal.valve(:,6));

ratioD.all = sitesY_valve;
ratioD.valve = ratioD.all(dorsal.valve(:,2)==1);
ratioD.novalve = ratioD.all(~(dorsal.valve(:,2)==1));

% vHPC ratio: -160 (deep) --> 1, -20 (superficial) --> 0
sitesY_valve = 1-double(ventral.shock(:,6));% Invert so that deep = 1

ratioV.all = sitesY_valve;  % deep = 1
ratioV.all = 1 - ratioV.all;
ratioV.valve = ratioV.all(ventral.valve(:,2)==1);
ratioV.novalve = ratioV.all(~(ventral.valve(:,2)==1));


x = [ones(size(ratioV.novalve)) ; ones(size(ratioV.valve))*2];
y = [ratioV.novalve ; ratioV.valve];
figure,
cdfplot(ratioV.novalve),hold on
cdfplot(ratioV.valve),hold on
% boxplot(y,x)
[h p] = ranksum(ratioV.novalve , ratioV.valve)


x = [ones(size(ratioD.novalve)) ; ones(size(ratioD.valve))*2];
y = [ratioD.novalve ; ratioD.valve];
figure,
cdfplot(ratioD.novalve),hold on
cdfplot(ratioD.valve),hold on
% boxplot(y,x)
[h p] = ranksum(ratioD.novalve , ratioD.valve)
end