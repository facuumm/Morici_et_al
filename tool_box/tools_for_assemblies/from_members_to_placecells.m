function from_members_to_placecells(path)

% HACER UN SCRIPT PARA PODER DETERMINAR QUE PC SON MIEMBROS DE QUE.
% LA IDEA ES TRATAR DE DEFINIR SI ESTAN TUNEADAS AL ESPACIO ALGUNOS MAS QUE
% OTROS.
% 
% 
% ----- Spike shape analysis routine -----
% Analyze one by one folder inside path.
% Lift spks and clusters to assign waveforms to clusters.
% Search for the waveform and calculates the spike duration
% of the spike to the middle of the valley (A) and the time between
% the valley and the second maximum (B). Save the figures and the values of A and B.
%
% syntax: waveform_parameters(path,sf)
%
% INPUTS
% path: Str, contains the path of the sessions you want to analyze.
%       This scripts get inside each subfolder of the introduced path and 
%       will get into the 'Spikesorting' folter where it will look for the
%       output of loadSpikes.m
% sf: int, sampling frequency. default 20kHz
%
% OUTPUS
% It will save the data in a mat file called 'Cell_Classification_using_waveform'
% It contains a matrix called 'A_B_waveforms' that contains data organaized as follows:
% Columns: 1st: Cluster_id / 2nd: cluster_channel / 3rd: A_distance / 4th: B_distance
% Also, It will save the plot aindividual waveforms with the parameters.
% Daniela Piña Novo 01/2018
% Camila Zold 09/19
% Facundo Morici 04/2023

sf = 20000; % sampling frequency
Segmento_BaseLineIndices = (1:1:10); % to calculate baseline (10)

%List of folders from the path
files = dir(path);
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags);
clear files dirFlags

for t = 1 : length(subFolders)-2
    session = [subFolders(t+2).folder,'\',subFolders(t+2).name,'\Spikesorting'];
    cd(session)
    
    
    x = dir([cd,'\dorsalventral_assemblies_reward*.mat']);
    load(x(end).name)
    clear x
    TH.reward = Th;
    clear Th criteria_fr criteria_n pat
    
    x = dir([cd,'\dorsalventral_assemblies_aversive*.mat']);
    load(x(end).name)
    clear x
    TH.aversive = Th;
    clear Th criteria_fr criteria_n pat
    
    % encocntrar una forma de definir cuales son dorsales y ventrales.
    
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
end

end
