% ----- Spike shape analysis routine -----
% Analyze one by one folder inside path.
% Lift spks and clusters to assign waveforms to clusters.
% Search for the waveform and calculates the spike duration
% of the spike to the middle of the valley (A) and the time between
% the valley and the second maximum (B). Save the figures and the values of A and B.

% This scripts get inside each subfolder of the introduced path and  will get into the
% 'Spikesorting' folter where it will look for the output of loadSpikes.m
% Output: '\Cell_Classification_using_waveform.mat'

% It contains a matrix called 'A_B_waveforms' that contains data organaized as follows:
% Columns: 1st: Cluster_id / 2nd: cluster_channel / 3rd: A_distance / 4th: B_distance

% Camila Zold 09/19
% Facundo Morici 04/2023



clear;
clc;
close all


%Parametros

path = {'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat103\usable';'E:\Rat132\recordings\in_pyr';'E:\Rat165\in_pyr'};

sf = 20000;
Segmento_BaseLineIndices = (1:1:10); %este es el n�mero de puntos que usa para calcular el baseline (10). cambiar de ser necesario

function wave_form_parameters(path)
    %List of folders from the path
    files = dir(path);
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
    clear files dirFlags
    
    for t = 1 : length(subFolders)-2
        A_B_waveforms = [];
        session = [subFolders(t+2).folder,'\',subFolders(t+2).name,'\Spikesorting'];
        cd(session)
        
        x = dir([cd,'\*spikes.cellinfo.mat']);
        load(x.name,'spikes')
        
        if ~exist('waveforms_figures', 'dir') %creation of folder to store waveforms
            mkdir waveforms_figures
        else % if it exist, delete files inside to overwrite
            delete('waveforms_figures/*.*')
        end
        
        for i = 1 : size(spikes.filtWaveform,2)
            
            waveform = spikes.filtWaveform{i};
            % polarity changment of inverted neurons
            x = and(and(min(waveform) < 0 , max(waveform) >0) , abs(min(waveform)) > max(waveform));
            if x
                waveform = waveform; %stay the same
            else
                waveform = waveform*-1; %change sign
            end
            clear x
            
            %plot spike
            figure(1); clf; hold on
            plot(waveform)
            xlabel('Samples','fontsize',10,'fontweight','bold')
            ylabel('Amplitude','fontsize',10,'fontweight','bold')
            hold on
            scatter([1:length(waveform)],waveform,'*r'); %todos los ptos
            pause (1);
            
            
            diff_waveform = diff(waveform); %first derived
            diff_waveform_round = round(diff_waveform);
            
            [Yvalue_min,Min_Maxchannel]=min(waveform); %minimo sobre el waveform
            
            %detection of sign changement to get max1 y max2
            Sign=sign(diff_waveform_round);
            
            
            
            for iind=2:Min_Maxchannel-4
                if (Sign(iind) ~= Sign(iind-1))
                    if ~exist('PtosCambioDeSigno1','var')
                        PtosCambioDeSigno1= iind;
                    else
                        PtosCambioDeSigno1= cat(1,PtosCambioDeSigno1,iind);
                    end
                end
            end
            clear iind
            
            for iind=(Min_Maxchannel+4):length(Sign)
                if (Sign(iind) ~= Sign(iind-1))
                    if ~exist('PtosCambioDeSigno2','var')
                        PtosCambioDeSigno2= iind;
                    else
                        PtosCambioDeSigno2= cat(1,PtosCambioDeSigno2,iind);
                    end
                end
            end
            clear iind
            
            if (exist('PtosCambioDeSigno1','var') && exist('PtosCambioDeSigno1','var'))
                Max1_Maxchannel=PtosCambioDeSigno1(end); %maximo 1 sobre la 1ra derivada
                Yvalue_max1=waveform(Max1_Maxchannel);
                clear PtosCambioDeSigno1
                
                Max2_Maxchannel=PtosCambioDeSigno2(1); %maximo 2 sobre la 1ra derivada
                Yvalue_max2=waveform(Max2_Maxchannel);
                clear PtosCambioDeSigno2
                
                %plot de los 3 puntos
                scatter(Max1_Maxchannel,Yvalue_max1,'o','g','filled','LineWidth',0.75) %maximo1
                scatter(Min_Maxchannel,Yvalue_min,'o','g','filled','LineWidth',0.75) %minimo
                scatter(Max2_Maxchannel,Yvalue_max2,'o','g','filled','LineWidth',0.75) %maximo2
                pause (1)
                
                %calcula distancia A=duraci�n del valle de la espiga en el 50% del valle; B=desde el valle (min) hasta el segundo m�ximo
                Distancia_B=Max2_Maxchannel-Min_Maxchannel;
                DistanciaB_ms=(Distancia_B/sf)*1000;
                %plot de la distancia B en lineas
                P1=[Min_Maxchannel Yvalue_min];P2=[Max2_Maxchannel Yvalue_min]; %linea horizontal
                P3=[Max2_Maxchannel Yvalue_min];P4=[Max2_Maxchannel Yvalue_max2]; %linea vertical
                plot([P1(1) P2(1)],[P1(2) P2(2)],'b','LineWidth',2) %linea horizontal--distancia C
                plot([P3(1) P4(1)],[P3(2) P4(2)],'b','linestyle', '--') %linea vertical--distancia C
                pause (1)
                
                
                % calcula BaseLine
                %Para el Baseline construyo un segmento tomando el numero de puntos que estan en parametros
                MediaSegmento=mean(waveform(Segmento_BaseLineIndices)); %Media en Y linea de base
                %plot Baseline
                P7=[length(waveform)-(length(waveform)-1) MediaSegmento];P8=[length(waveform) MediaSegmento]; %linea horizontal--baseline
                plot([P7(1) P8(1)],[P7(2) P8(2)],'r','LineWidth',1.2) %linea horizontal baseline
                pause (1)
                
                %calcula distancia A
                MediaAltura= (MediaSegmento+Yvalue_min)/2;
                Pto1_ancho_antes=find(waveform(Max1_Maxchannel:Min_Maxchannel)>MediaAltura,1,'last');
                Pto1_ancho_dspues=find(waveform(Max1_Maxchannel:Min_Maxchannel)<MediaAltura,1,'first');
                Pto1_ancho_antes_correg=(Max1_Maxchannel+Pto1_ancho_antes)-1;
                Pto1_ancho_dspues_correg=(Max1_Maxchannel+Pto1_ancho_dspues)-1;
                
                Pto2_ancho_antes=find(waveform(Min_Maxchannel:Max2_Maxchannel)<MediaAltura,1,'last');
                Pto2_ancho_dspues=find(waveform(Min_Maxchannel:Max2_Maxchannel)>MediaAltura,1,'first');
                Pto2_ancho_antes_correg=(Min_Maxchannel+Pto2_ancho_antes)-1;
                Pto2_ancho_dspues_correg=(Min_Maxchannel+Pto2_ancho_dspues)-1;
                
                
                %interpolo primer pto para calcular ancho
                Y0=waveform(Pto1_ancho_antes_correg);
                Y1=waveform(Pto1_ancho_dspues_correg);
                Y=MediaAltura;
                
                X0=Pto1_ancho_antes_correg;
                X1=Pto1_ancho_dspues_correg;
                
                m=(Y1-Y0)/(X1-X0);
                n=Y0-(m*X0);
                X_primerPto=abs((Y-n)/m);
                clear Y1 Y0 Y X1 X0 n m
                
                %interpolo segundo pto para calcular ancho
                Y0=waveform(Pto2_ancho_antes_correg);
                Y1=waveform(Pto2_ancho_dspues_correg);
                Y=MediaAltura;
                
                X0=Pto2_ancho_antes_correg;
                X1=Pto2_ancho_dspues_correg;
                
                m=(Y1-Y0)/(X1-X0);
                n=Y0-(m*X0);
                X_segundoPto=abs((Y-n)/m);
                clear Y1 Y0 Y X1 X0 n m
                
                Distancia_A=X_segundoPto-X_primerPto;
                DistanciaA_ms=(Distancia_A/sf)*1000;
                
                %plot de la distancia A en lineas
                P5=[X_primerPto MediaAltura];P6=[X_segundoPto MediaAltura]; %linea horizontal distancia D
                plot([P5(1) P6(1)],[P5(2) P6(2)],'b','LineWidth',2) %distancia D
                %     title({['File:', ];['Cluster:', num2str(SpikesToLoad(clust))]; ['Distancia A (UM): ',num2str(Distancia_A),'(',num2str(DistanciaA_ms),'ms)'];['Distancia B (UM): ',num2str(Distancia_B),'(',num2str(DistanciaB_ms),'ms)']},'fontsize',10); hold on
                pause (1.5)
                
                details = [spikes.cluID(i) , spikes.maxWaveformCh1(i)];
                A_B_waveforms=[A_B_waveforms ; details , DistanciaA_ms,DistanciaB_ms];
                clear Max1_Maxchannel Min_Maxchannel
%             else
%                 Distancia_A = 333;
%                 DistanciaA_ms = 333;
%                 Distancia_B = 333;
%                 DistanciaB_ms = 333;
%                 
%                 details = [spikes.cluID(i) , spikes.maxWaveformCh1(i)];
%                 A_B_waveforms = [A_B_waveforms ; details , DistanciaA_ms,DistanciaB_ms];
%                 clear Max1_Maxchannel Min_Maxchannel
            end
            FigName=['Cluster# ',num2str(details(1)), '.fig'];
            savefig([cd,'\waveforms_figures\',FigName])
            clear details
        end
        save([cd,'\Cell_Classification_using_waveform.mat'],'A_B_waveforms')
        clear A_B_waveforms
        t
    end
    tt
end
    




%% From here 
todo = [];
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
        load('Cell_Classification_using_waveform.mat')
        load('burst_index.mat')
        
        if isequal(A_B_waveforms(:,1) , burst_index(:,1))
            todo = [todo ; A_B_waveforms , burst_index(:,2)];
        else
            break
            t
            tt
        end
        clear A_B_waveforms burst_index
    end
end
    


[idx,C] = kmeans(todo(:,3:4),3,'Start','cluster');

int = idx==2;
pyr = idx==3;
dHPC = todo(:,2)>63;
vHPC = todo(:,2)<=63;

subplot(131)
scatter3(todo(pyr,3),todo(pyr,4),todo(pyr,5),'filled','MarkerFaceColor',[0.6 0.0 0.6]),hold on%,'MarkerFaceAlpha',0.5),hold on
scatter3(todo(int,3),todo(int,4),todo(int,5),'filled','MarkerFaceColor',[0.1986 0.7214 0.6310])%,'MarkerFaceAlpha',0.5)
ylim([0 1]), xlim([0 0.5]), zlim([0 8])
xlabel('Distance A')
ylabel('Distance B')
zlabel('Burst Index')
title('Pyr/Int classification - All')

% dHPC
pyr_dHPC = and(pyr,dHPC);
int_dHPC = and(int,dHPC);


%vHPC
pyr_vHPC = and(pyr,vHPC);
int_vHPC = and(int,vHPC);


subplot(132)
scatter3(todo(pyr_dHPC,3),todo(pyr_dHPC,4),todo(pyr_dHPC,5),'filled','MarkerFaceColor',[0.6 0.0 0.6]),hold on%,'MarkerFaceAlpha',0.5),hold on
scatter3(todo(int_dHPC,3),todo(int_dHPC,4),todo(int_dHPC,5),'filled','MarkerFaceColor',[0.1986 0.7214 0.6310])%,'MarkerFaceAlpha',0.5)
ylim([0 1]), xlim([0 0.5]), zlim([0 8])
xlabel('Distance A')
ylabel('Distance B')
zlabel('Burst Index')
title('Pyr/Int classification - dHPC')

subplot(133)
scatter3(todo(pyr_vHPC,3),todo(pyr_vHPC,4),todo(pyr_vHPC,5),'filled','MarkerFaceColor',[0.6 0.0 0.6]),hold on%,'MarkerFaceAlpha',0.5),hold on
scatter3(todo(int_vHPC,3),todo(int_vHPC,4),todo(int_vHPC,5),'filled','MarkerFaceColor',[0.1986 0.7214 0.6310])%,'MarkerFaceAlpha',0.5)
ylim([0 1]), xlim([0 0.5]), zlim([0 8])
xlabel('Distance A')
ylabel('Distance B')
zlabel('Burst Index')
title('Pyr/Int classification - vHPC')



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
        load('Cell_Classification_using_waveform.mat')
        load('burst_index.mat')
        
        for i = 1:length(burst_index)
            if i ==1
            cond1 = and(todo(:,1) == A_B_waveforms(i,1) , todo(:,2) == A_B_waveforms(i,2));
            cond2 = and(todo(:,3) == A_B_waveforms(i,3) , todo(:,4) == A_B_waveforms(i,4));
            else
            cond1 = cond1 + and(todo(:,1) == A_B_waveforms(i,1) , todo(:,2) == A_B_waveforms(i,2));
            cond2 = cond2 + and(todo(:,3) == A_B_waveforms(i,3) , todo(:,4) == A_B_waveforms(i,4));
            end
        end

        tabulador= and(logical(cond1) , logical(cond2));
        clear cond1 cond2 cond3

        Cell_type_classification = [todo(tabulador,:) , pyr(tabulador,:) , int(tabulador,:)];
        save('Cell_type_classification.mat' , 'Cell_type_classification')
        clear Cell_type_classification
    end
end

