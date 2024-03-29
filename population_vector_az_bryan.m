%% Script to calculate population vector correlation of 2D ratemaps 
% Sessions need to have at least 5 neurons to be included in the pv corr
% analysis
% Silva 2023

%%
clear
clc
close all
%% Parameters
path = {'\\Maryjackson\e\Rat127\Ephys\pyr';'\\Maryjackson\e\Rat128\Ephys\in_pyr\ready';'\\Maryjackson\e\Rat103\usable';'\\Maryjackson\e\Rat132\recordings\in_pyr'; '\\Maryjackson\e\Rat165\in_pyr'};%List of folders from the path
% path = {'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat103\usable';'E:\Rat132\recordings\in_pyr'};%List of folders from the path
%% Output 
pv_corr.dhpc = []; 
pv_corr.vhpc = []; 
%%
doplot=1
%% Accumulated population vector correlation 
for tt = 1:length(path)
    
    %List of folders from the path
    files = dir(path{tt});
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
    clear files dirFlags
    
    temp_dhpc= [];
    temp_vhpc= [];
    
    for t = 1 : length(subFolders)-2
        disp(['-- Initiating analysis of folder #' , num2str(t) , ' from rat #',num2str(tt) , ' --'])
        session = [subFolders(t+2).folder,'\',subFolders(t+2).name];
        cd(session)
        cd 'Spikesorting'
        
        temp = [0 0];
        if isfile('dHPC_pc.mat')
            load('dHPC_pc.mat');
            temp(1) = size(dHPC,2);
        end
        if isfile('vHPC_pc.mat')
            load('vHPC_pc.mat');
            temp(2) = size(vHPC,2);
        end
        min_neurons = min(temp);
        
        if min_neurons==0 | min_neurons <5
            continue
        end
        pv_ave_dhpc=[];
        pv_rew_dhpc=[];
        if isfile('dHPC_pc.mat')
            if size(dHPC,2)>=5 
                %Create rate map stacks 
                for n=1:min_neurons%size(dHPC,2)
                    pv_ave_dhpc(:,:,n) = dHPC{n}.frMap_ave; 
                    pv_rew_dhpc(:,:,n) = dHPC{n}.frMap_rew;
                end
                %%
                % Vector bin correlation - dorsal
                sesion_corr = NaN(60,1);
                for c=1:size(pv_ave_dhpc,2) %Iterate through all the bins
                    %Create a bin vector (with the data of all neurons in that bin)
                    pv_ave = squeeze(pv_ave_dhpc(1,c,:));
            
                    pv_rew = squeeze(pv_rew_dhpc(1,c,:));
            
                    %Correlation between ave vector bin and rew vector bin
                    cor = corrcoef(pv_ave, pv_rew,'Rows','pairwise');
                    sesion_corr(c,1) = cor(1,2);    
                end 
                temp_dhpc = [temp_dhpc;sesion_corr];
                
%                 clear pv_ave_dhpc pv_rew_dhpc pv_ave pv_rew
                
                
            end
        end 
        
        pv_ave_vhpc=[];
        pv_rew_vhpc=[];
        if isfile('vHPC_pc.mat')
            load('vHPC_pc.mat');
            if size(vHPC,2)>=5
                %Create rate map stacks 
                for n=1:min_neurons%size(vHPC,2)
                    pv_ave_vhpc(:,:,n) = vHPC{n}.frMap_ave; 
                    pv_rew_vhpc(:,:,n) = vHPC{n}.frMap_rew;
                end 
                
                % Vector bin correlation - ventral
                sesion_corr = NaN(60,1);
                for c=1:size(pv_ave_vhpc,2) %Iterate through all the bins
                    %Create a bin vector (with the data of all neurons in that bin)
                    pv_ave = squeeze(pv_ave_vhpc(1,c,:));
            
                    pv_rew = squeeze(pv_rew_vhpc(1,c,:));
                    
                    %Correlation between ave vector bin and rew vector bin
                    cor = corrcoef(pv_ave, pv_rew,'Rows','pairwise');
                    sesion_corr(c,1) = cor(1,2);    
                end 
                
                
                if doplot==1
                    figure(1)
                    subplot(2,2,1)
                    imagesc(squeeze(pv_ave_dhpc)')
                    ylabel('Neurons')
                    subplot(2,2,2)
                    imagesc(squeeze(pv_rew_dhpc)')
                    xlabel('Space')
                    figure(1)
                    subplot(2,2,3)
                    imagesc(squeeze(pv_ave_vhpc)')
                    ylabel('Neurons')
                    subplot(2,2,4)
                    imagesc(squeeze(pv_rew_vhpc)')
                    xlabel('Space')
                    pause()
                end
                
                temp_vhpc = [temp_vhpc;sesion_corr];  
                
                clear pv_ave_vhpc pv_rew_vhpc pv_ave pv_rew
            end
        end 
              
    end
    
    pv_corr.dhpc = [pv_corr.dhpc;temp_dhpc];
    pv_corr.vhpc = [pv_corr.vhpc;temp_vhpc];
    disp(['-------- Finished rat#' , num2str(tt) , ' --------'])
    disp(' ')
end %iterate rats 
% Plots: 
figure(1);clf;hold on; 
title('Population vector correlation')
h = cdfplot(pv_corr.dhpc);hold on;
h.Color = 'green';
h=cdfplot(pv_corr.vhpc);hold on;
h.Color = 'blue';
ylabel('Cumulative freq');
xlabel('Spatial bin corr. coef.(r)')


%% Population vector per session 
figure(2);clf;hold on; 
for tt = 1:length(path)
    %List of folders from the path
    files = dir(path{tt});
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
    clear files dirFlags
    
    temp_dhpc= [];
    temp_vhpc= [];
    
    for t = 1 : length(subFolders)-2
        disp(['-- Initiating analysis of folder #' , num2str(t) , ' from rat #',num2str(tt) , ' --'])
        session = [subFolders(t+2).folder,'\',subFolders(t+2).name];
        cd(session)
        cd 'Spikesorting'
        
        if isfile('dHPC_pc.mat')
            load('dHPC_pc.mat');
            if size(dHPC,2)>=10 
                %Create rate map stacks 
                for n=1:size(dHPC,2)
                    pv_ave_dhpc(:,:,n) = dHPC{n}.frMap_ave; 
                    pv_rew_dhpc(:,:,n) = dHPC{n}.frMap_rew;
                end
             
                % Vector bin correlation - dorsal
                sesion_corr = NaN(60,1);
                for c=1:size(pv_ave_dhpc,2) %Iterate through all the bins
                    %Create a bin vector (with the data of all neurons in that bin)
                    pv_ave = squeeze(pv_ave_dhpc(1,c,:));
            
                    pv_rew = squeeze(pv_rew_dhpc(1,c,:));
            
                    %Correlation between ave vector bin and rew vector bin
                    cor = corrcoef(pv_ave, pv_rew);
                    sesion_corr(c,1) = cor(1,2);  
                end 
                h = cdfplot(sesion_corr);hold on;
                h.Color = 'green';
                clear pv_ave_dhpc pv_rew_dhpc pv_ave pv_rew
                
            end
        end 
        
        if isfile('vHPC_pc.mat')
            load('vHPC_pc.mat');
            if size(vHPC,2)>=10
                %Create rate map stacks 
                for n=1:size(vHPC,2)
                    pv_ave_vhpc(:,:,n) = vHPC{n}.frMap_ave; 
                    pv_rew_vhpc(:,:,n) = vHPC{n}.frMap_rew;
                end 
                
                % Vector bin correlation - dorsal
                sesion_corr = NaN(60,1);
                for c=1:size(pv_ave_vhpc,2) %Iterate through all the bins
                    %Create a bin vector (with the data of all neurons in that bin)
                    pv_ave = squeeze(pv_ave_vhpc(1,c,:));
            
                    pv_rew = squeeze(pv_rew_vhpc(1,c,:));
            
                    %Correlation between ave vector bin and rew vector bin
                    cor = corrcoef(pv_ave, pv_rew);
                    sesion_corr(c,1) = cor(1,2);    
                    
                end 
                h = cdfplot(sesion_corr);hold on;
                h.Color = 'blue';   
                clear pv_ave_vhpc pv_rew_vhpc pv_ave pv_rew
            end
        end 
              
    end
    
    disp(['-------- Finished rat#' , num2str(tt) , ' --------'])
    disp(' ')
end %iterate rats 
ylabel('Cumulative freq');
xlabel('Spatial bin corr. coef.(r)')
title('Population vector correlation')