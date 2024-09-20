%% General variables
criteria_fr = 0;
criteria_type = 0; % criteria for celltype (0:pyr, 1:int, 2:all)

%% Load LFP
load('lfp.mat')

%% Loading TS of the sessions
disp('Uploading session time stamps')
load('session_organization.mat')
load('behavioral_data.mat')

baselineTS = baselineTS./1000;                  aversiveTS = aversiveTS./1000;                    rewardTS = rewardTS./1000;
aversiveTS_run = aversiveTS_run./1000;          rewardTS_run = rewardTS_run./1000;

%% Load ripples
if exist('ripplesD_customized2.csv')
    ripplesD = table2array(readtable('ripplesD_customized2.csv'));
    RD = true;
else
    RD = false;
end

if exist('ripplesV_customized2.csv')
    ripplesV = table2array(readtable('ripplesV_customized2.csv'));
    RV = true;
else
    RV = false;
end

if and(RD,RV)
    RB = true;
    % coordination
    coordinated = [];
    coordinatedV = [];
    cooridnated_event = [];
    for i = 1:length(ripplesD)
        r = ripplesD(i,:);
        tmp = sum(and(ripplesV(:,2)>= r(1,2)-0.1, ripplesV(:,2)<= r(1,2)+0.1));
        if tmp>0
            z = ripplesV(and(ripplesV(:,2)>= r(1,2)-0.1, ripplesV(:,2)<= r(1,2)+0.1),:);
            coordinatedV = [coordinatedV ; z];
            [p,indice] = min(abs(r(2)-z(:,2)));
            coordinated = [coordinated ; r];
            
            peak = min(min(z(indice,2),r(2)))+abs(z(indice,2)-r(2))/2;
            low = min([r(1) , z(indice,1)]);
            up = max([r(3) , z(indice,3)]);
            cooridnated_event = [cooridnated_event ; low , peak , up];
            
            clear tmp2 tmp1 p indice z peak low up
        end
        clear r
    end
    clear x tmp i
    
    [C,IA,IC] = unique(coordinatedV(:,1));
    coordinatedV  = coordinatedV(IA,:); clear C IA IC
    % Store events time stamps
    ripple_event.all = Restrict(cooridnated_event,aversiveTS);
end

%% load spks
disp('Uploading Spiking activity')
cd 'Spikesorting'
[clusters , numberD , numberV , spks , spks_dHPC , spks_vHPC , cellulartype] = load_SU_FM(cd,criteria_type,criteria_fr,aversiveTS_run,rewardTS_run);


for i = 10 : 30 %length(ripple_event.all)
    idx = [ripple_event.all(i,2)-0.2 ripple_event.all(i,2)+0.2];
    
    % Restrictitng lfp
    lfp1 = Restrict(dHPC,idx);
    lfp2 = Restrict(vHPC1,idx);
    
    figure
    subplot(311)
    plot(lfp1(:,1) , lfp1(:,2)),hold on
    plot(lfp2(:,1) , lfp2(:,2)),hold on
    
    % filtering
    lfp1 = FilterLFP(lfp1,'passband',[100 250]);
    lfp2 = FilterLFP(lfp2,'passband',[100 250]);
    subplot(312)
    plot(lfp1(:,1) , lfp1(:,2)),hold on
    plot(lfp2(:,1) , lfp2(:,2)),hold on    

    % SUs
    subplot(313)
    iii=0
for ii = 1 : length(clusters.all)
    cluster = clusters.all(ii,1);
%     celltype = logical(cellulartype(cellulartype(:,1) == cluster,2));
%     if celltype
        iii = iii+1;
        s = spks(spks(:,1)==cluster,2);
        s = Restrict(s,idx);
        s = s';
        xspikes = repmat(s,3,1);
        yspikes = nan(size(xspikes));
        
        if not(isempty(yspikes))
            yspikes(1,:) = iii-1;
            yspikes(2,:) = iii;
        end
        
        %                             if i >= C
        if ii>length(clusters.vHPC)
            
            plot(xspikes,yspikes,'Color','r','LineWidth',1),hold on
        else
            plot(xspikes,yspikes,'Color','k','LineWidth',1),hold on
        end
        %                             yline(ii)
%     end
end
xlim(idx)
ylim([0 length(clusters.all)])
ylabel('Neuron#')    
    
    
    
end
