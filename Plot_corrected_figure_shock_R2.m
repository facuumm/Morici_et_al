%% dHPS SUs
tmp1 = [];
response.pre.dHPC = [];
for i = 1 : size(Pre.reward.dHPC,2)
    [~ , ii] = min(abs(T-(-0.20)));    [~ , iii] = min(abs(T-(0.20)));
    
    t = Pre.reward.dHPC(:,i);
    
    M = nanmean([t(1:ii) ; t(iii:end)]);
    S = nanstd([t(1:ii) ; t(iii:end)]);
    
    t = (t./M);
    t = Smooth(t,2);
    
    [~ , ii] = min(abs(T-(-0.15)));    [~ , iii] = min(abs(T-(0.15)));
    response.pre.dHPC = [response.pre.dHPC ; max(t(ii:iii))];
    
    tmp1 = [tmp1 , t]; clear M S
    
end

tmp2 = [];
response.post.dHPC = [];
for i = 1 : size(Post.reward.dHPC,2)
    [~ , ii] = min(abs(T-(-0.20)));    [~ , iii] = min(abs(T-(0.20)));
    
    t = Post.reward.dHPC(:,i);
    
    M = nanmean([t(1:ii) ; t(iii:end)]);
    S = nanstd([t(1:ii) ; t(iii:end)]);
    
    t = (t./M);
    t = Smooth(t,2);
    
    [~ , ii] = min(abs(T-(-0.15)));    [~ , iii] = min(abs(T-(0.15)));
    response.post.dHPC = [response.post.dHPC ; max(t(ii:iii))];
    
    tmp2 = [tmp2 , t]; clear M S
    
end



figure,
i = and(I.dHPC(:,1) == 1 , I.dHPC(:,5)==1);
% i = I.dHPC(:,1) == 1;
subplot(221),
m = nanmean((tmp1(:,i))');
s = nansem((tmp1(:,i))');
plot(T,m,'k'),hold on
ciplot(m-s , m+s , T , 'k'),alpha 0.5

m = nanmean((tmp2(:,i))');
s = nansem((tmp2(:,i))');
plot(T,m,'r'),hold on
ciplot(m-s , m+s , T , 'r'),alpha 0.5
xlim([-0.3 0.3])
ylim([0 6])

[h p] = signrank(response.pre.dHPC(i) , response.post.dHPC(i),'tail','left')


subplot(222),

% Downsampling the data of non-shock responssivness
x = [1:length(i)];
i = and(not(I.dHPC(:,1) == 1) , I.dHPC(:,5)==1);
x = x((i));
ii = randperm(length(x));
x = x(ii(1:sum(i)));

m = nanmean((tmp1(:,x))');
s = nansem((tmp1(:,(i)))');
plot(T,m,'k'),hold on
ciplot(m-s , m+s , T , 'k'),alpha 0.5

m = nanmean((tmp2(:,x))');
s = nansem((tmp2(:,x))');
plot(T,m','r'),hold on
ciplot(m-s , m+s , T , 'r'),alpha 0.5
xlim([-0.3 0.3])
ylim([0 6])

[h p] = signrank(response.pre.dHPC(i) , response.post.dHPC(i),'tail','left')

%% vHPS SUs
tmp1 = [];
response.pre.vHPC = [];
for i = 1 : size(Pre.aversive.vHPC,2)
    [~ , ii] = min(abs(T-(-0.2)));    [~ , iii] = min(abs(T-(0.2)));
    
    t = Pre.aversive.vHPC(:,i);
    
    M = nanmean([t(1:ii) ; t(iii:end)]);
    S = nanstd([t(1:ii) ; t(iii:end)]);
    
    t = (t./M);
    t = Smooth(t,2);
    
    [~ , ii] = min(abs(T-(-0.15)));    [~ , iii] = min(abs(T-(0.15)));
    response.pre.vHPC = [response.pre.vHPC ; max(t(ii:iii))];
    
    tmp1 = [tmp1 , t]; clear M S
    
end

tmp2 = [];
response.post.vHPC = [];
for i = 1 : size(Post.aversive.vHPC,2)
    [~ , ii] = min(abs(T-(-0.2)));    [~ , iii] = min(abs(T-(0.2)));
    
    t = Post.aversive.vHPC(:,i);
    
    M = nanmean([t(1:ii) ; t(iii:end)]);
    S = nanstd([t(1:ii) ; t(iii:end)]);
    
    t = (t./M);
    t = Smooth(t,2);
    
    [~ , ii] = min(abs(T-(-0.15)));    [~ , iii] = min(abs(T-(0.15)));
    response.post.vHPC = [response.post.vHPC ; max(t(ii:iii))];
    
    tmp2 = [tmp2 , t]; clear M S
    
end

i = and(I.vHPC(:,1) == 1 , I.vHPC(:,5)==1);
subplot(223),
m = nanmean((tmp1(:,i))');
s = nansem((tmp1(:,i))');
plot(T,m,'k'),hold on
ciplot(m-s , m+s , T , 'k'),alpha 0.5

m = nanmean((tmp2(:,i))');
s = nansem((tmp2(:,i))');
plot(T,m,'r'),hold on
ciplot(m-s , m+s , T , 'r'),alpha 0.5
xlim([-0.3 0.3])
ylim([0 6])

[h p] = signrank(response.pre.vHPC(i) , response.post.vHPC(i),'tail','left')

% Downsampling the data of non-shock responssivness
x = [1:length(i)];
i = and(not(I.vHPC(:,1) == 1) , I.vHPC(:,5)==1);
x = x((i));
ii = randperm(length(x));
x = x(ii(1:sum(i)));

subplot(224),
m = nanmean((tmp1(:,x))');
s = nansem((tmp1(:,x))');
plot(T,m,'k'),hold on
ciplot(m-s , m+s , T , 'k'),alpha 0.5

m = nanmean((tmp2(:,x))');
s = nansem((tmp2(:,x))');
plot(T,m','r'),hold on
ciplot(m-s , m+s , T , 'r'),alpha 0.5
xlim([-0.3 0.3])
ylim([0 6])
[h p] = signrank(response.pre.vHPC(i) , response.post.vHPC(i),'tail','left')



%% Post-Pre vs Shock correlation
response.delta.dHPC = response.post.dHPC - response.pre.dHPC;
response.delta.vHPC = response.post.vHPC - response.pre.vHPC;

x = response.shock.dHPC(and((I.dHPC(:,1) == 1) , I.dHPC(:,5)==1));
y = response.delta.dHPC(and((I.dHPC(:,1) == 1) , I.dHPC(:,5)==1));
figure,fitlm(x,y)
subplot(121),scatter(x,y,'filled'), ylim([-1.5 2]) , xlim([-0.5 5])
subplot(122),plot(ans), ylim([-1.5 2]) , xlim([-0.5 5])

x = response.shock.vHPC(and((I.vHPC(:,1) == 1) , I.vHPC(:,5)==1));
y = response.delta.vHPC(and((I.vHPC(:,1) == 1) , I.vHPC(:,5)==1));
x(43) = []; x(18) = []; y(43) = []; y(18) = []; 
figure,fitlm(x,y)
subplot(121),scatter(x,y,'filled'), ylim([-1.5 2]) , xlim([-0.5 5])
subplot(122),plot(ans), ylim([-1.5 2]) , xlim([-0.5 5])

%% Correlation Reactivation vs Gain
C = [];
for i = 1 : size(SAVED.aversive,2)
    tmp = SAVED.aversive{i};
    tmp1 = [I.vHPC(:,4) I.vHPC(:,3)];
    tmp1 = and(tmp1(:,1) == tmp.Session(1) , tmp1(:,2) == tmp.Session(2));
    TMP1 = and((I.vHPC(:,1) == 1) , I.vHPC(:,5)==1); tmp1 = and(TMP1,tmp1);
    tmp1 = I.vHPC(tmp1,2);
    
    tmpX = [I.vHPC(:,4) I.vHPC(:,3)];
    tmp2 = response.post.vHPC(and(tmpX(:,1) == tmp.Session(1) , tmpX(:,2) == tmp.Session(2)),1);
%     tmp3 = response.pre.vHPC(and(tmpX(:,1) == tmp.Session(1) , tmpX(:,2) == tmp.Session(2)),1);
%     
%     tmp2 = (tmp2-tmp3)./mean([tmp2,tmp3],2);
    
    for ii = 1 : size(tmp.ids,1)
        if ismember(tmp.ids(ii),tmp1)
            index = tmp1==tmp.ids(ii);
            C = [C  ; tmp.Reactivation tmp2(index)];
        end
    end
end

C(33,:) = [];
fitlm(C(:,1),C(:,2))
figure,plot(ans),ylim([0 10]),xlim([-0.6 1.2])
figure,scatter(C(:,1),C(:,2),'filled'),ylim([0 10]),xlim([-0.6 1.2])


C = [];
for i = 1 : size(SAVED.aversive,2)
    tmp = SAVED.aversive{i};
    tmp1 = [I.dHPC(:,4) I.dHPC(:,3)];
    tmp1 = and(tmp1(:,1) == tmp.Session(1) , tmp1(:,2) == tmp.Session(2));
    TMP1 = and((I.dHPC(:,1) == 1) , I.dHPC(:,5)==1); tmp1 = and(TMP1,tmp1);
    tmp1 = I.dHPC(tmp1,2);
    
    tmp2 = [I.dHPC(:,4) I.dHPC(:,3)];
    tmp2 = response.post.dHPC(and(tmp2(:,1) == tmp.Session(1) , tmp2(:,2) == tmp.Session(2)),1);
%     tmp3 = response.pre.dHPC(and(tmp2(:,1) == tmp.Session(1) , tmp2(:,2) == tmp.Session(2)),1);
%     
%     tmp2 = tmp2-tmp3;
    
    for ii = 1 : size(tmp.ids,1)
        if ismember(tmp.ids(ii),tmp1)
            index = tmp1==tmp.ids(ii);
            C = [C  ; tmp.Reactivation tmp2(index)];
        end
    end
end
fitlm(C(:,1),C(:,2))
figure,plot(ans),ylim([0 10]),xlim([-0.6 1.2])
figure,scatter(C(:,1),C(:,2),'filled'),ylim([0 10]),xlim([-0.6 1.2])

%% Comparison members vs no-members
xD = response.shock.dHPC(and((I.dHPC(:,1) == 1) , I.dHPC(:,5)==1));
yD = response.shock.dHPC(and((I.dHPC(:,1) == 1) , I.dHPC(:,5)==0));
xV = response.shock.vHPC(and((I.vHPC(:,1) == 1) , I.vHPC(:,5)==1));
yV = response.shock.vHPC(and((I.vHPC(:,1) == 1) , I.vHPC(:,5)==0));

% % subsampling no-members
% i = randperm(length(yD)); i = i(1:length(xD));
% yD = yD(i);
% i = randperm(length(yV)); i = i(1:length(xV));
% yV = yV(i);


y = [xD ; yD ; xV ; yV];
grps = [ones(size(xD,1),1) ; ones(size(yD,1),1)*2 ; ones(size(xV,1),1)*3 ; ones(size(yV,1),1)*4];

figure
scatter(grps,y,'filled','jitter','on', 'jitterAmount',0.1),hold on
scatter([1 2 3 4] , [nanmean(xD) ; nanmean(yD) ; nanmean(xV) ; nanmean(yV)],'filled'),xlim([0 5])

grps = [ones(size(xD,1),1) ; ones(size(yD,1),1) ; ones(size(xV,1),1)*2 ; ones(size(yV,1),1)*2];
grps = [grps,[ones(size(xD,1),1) ; ones(size(yD,1),1)*2 ; ones(size(xV,1),1) ; ones(size(yV,1),1)*2]];
[~,~,stats] = anovan(y,grps,"Model","interaction","Varnames",["structure","member"])
[results,~,~,gnames] = multcompare(stats,"Dimension",[1 2]);

%% percentage of shock neurons in joint assemblies

totalD = sum(I.dHPC(:,5)==1);
totalV = sum(I.vHPC(:,5)==1);
shockD = sum(and(I.dHPC(:,5)==1 , I.dHPC(:,1)==1));
shockV = sum(and(I.vHPC(:,5)==1 , I.vHPC(:,1)==1));
shockD = (shockD/totalD)*100; shockV = (shockV/totalV)*100;

figure
subplot(211),bar([shockD 100])
subplot(212),bar([shockV 100])

%% Decorrelation of speed and shock response
time = [-5:0.1:5];
[h hh] = min(abs(time-(1)));
decorrelated.dHPC = [];
for i = 1 : size(curve.dHPC,1)
    id = curve.id.dHPC(i,:);
    template = curve.speed.id;
    template = and(id(3)==template(:,1) , id(4)==template(:,2));
    [deco] = decorr(curve.dHPC(i,:)',curve.speed.data(:,template));
    n = nanmean(deco(1:hh));
    decorrelated.dHPC = [decorrelated.dHPC , deco]; clear deco n
end

decorrelated.vHPC = [];
for i = 1 : size(curve.vHPC,1)
    id = curve.id.vHPC(i,:);
    template = curve.speed.id;
    template = and(id(3)==template(:,1) , id(4)==template(:,2));
    [deco] = decorr(curve.vHPC(i,:)',curve.speed.data(:,template));
    n = nanmean(deco(1:hh));
    decorrelated.vHPC = [decorrelated.vHPC , deco]; clear deco n
end

%% Plot
% vHPC
% subplot(122)
i = curve.id.vHPC(:,2) == 1;
% m = nanmean(decorrelated.vHPC(:,i)');
% s = nansem(decorrelated.vHPC(:,i)');
% plot([-5:0.1:5],m,'g'), hold on
% ciplot(m-s , m+s , [-5:0.1:5] , 'g'),alpha 0.10
% xline(0),
% xline(1)

[~ , iii] = min(abs(time-0.1));
[~ , iiii] = min(abs(time-1));
values.vHPC.decorrelated = max(decorrelated.vHPC(iii:iiii,i))';

i = curve.id.vHPC(:,2) == 1;
x = curve.vHPC';% - nanmean(curve.vHPC(:,1:hh)');%normalization to baseline as R2 asked
% m = nanmean(x(i,:));
% s = nansem(x(i,:));
% plot([-5:0.1:5],m,'y'), hold on
% ciplot(m-s , m+s , [-5:0.1:5] , 'y'),alpha 0.1

values.vHPC.shock = [max(x(iii:iiii,i))'];
values.vHPC.all = [max(x(iii:iiii,:))'];
% 
% i = not(curve.id.vHPC(:,2) == 1);
% m = nanmean(x(i,:));
% s = nansem(x(i,:));
% plot([-5:0.1:5],m,'k'), hold on
% ciplot(m-s , m+s , [-5:0.1:5] , 'k'),alpha 0.1
% ylim([-0.25 1.1])
% 
% values.vHPC.no =  max(x(i,iii:iiii)')'; clear x



% dHPC
% subplot(121)
i = curve.id.dHPC(:,2) == 1;
% m = nanmean(decorrelated.dHPC(:,i)');
% s = nansem(decorrelated.dHPC(:,i)');
% plot(time,m,'g'), hold on
% ciplot(m-s , m+s , time , 'g'),alpha 0.10

values.dHPC.decorrelated =  max(decorrelated.dHPC(iii:iiii,i))';

i = curve.id.dHPC(:,2) == 1;
x = curve.dHPC';% - nanmean(curve.dHPC(:,1:hh)'); %normalization to baseline as R2 asked
% m = nanmean(x(i,:));
% s = nansem(x(i,:));
% plot(time,m,'y'), hold on
% ciplot(m-s , m+s , time , 'y'),alpha 0.1

values.dHPC.shock = max(x(iii:iiii,i))';
values.dHPC.all = max(x(iii:iiii,:))';

% i = not(curve.id.dHPC(:,2) == 1);
% m = nanmean(x(i,:));
% s = nansem(x(i,:));
% plot(time,m,'k'), hold on
% ciplot(m-s , m+s , time , 'k'),alpha 0.1
% ylim([-0.25 1.1])
% xline(0),
% xline(1)

% values.dHPC.no = max(x(i,iii:iiii)')'; clear c


%% Plot and Analyse
figure,
grps = [ones(size(values.dHPC.shock,1),1) ; ones(size(values.dHPC.decorrelated,1),1)*2 ; ones(size(values.vHPC.shock,1),1)*3 ; ones(size(values.vHPC.decorrelated,1),1)*4];
y = [values.dHPC.shock ; values.dHPC.decorrelated ; values.vHPC.shock ; values.vHPC.decorrelated]; 
scatter(grps,y,'filled','jitter','on', 'jitterAmount',0.1),ylim([-0.5 6]) , xlim([0 5]),hold on
scatter([1 2 3 4] , [nanmean(values.dHPC.shock) ; nanmean(values.dHPC.decorrelated) ; nanmean(values.vHPC.shock) ; nanmean(values.vHPC.decorrelated)] , 'filled')

% Run analysis
y = [values.dHPC.shock ; values.dHPC.decorrelated ; values.vHPC.shock ; values.vHPC.decorrelated]; 

structure = [ones(size(values.dHPC.shock,1),1) ; ones(size(values.dHPC.decorrelated,1),1) ; ...
             2*ones(size(values.vHPC.shock,1),1) ; 2*ones(size(values.vHPC.decorrelated,1),1)];

condition = [ones(size(values.dHPC.shock,1),1) ; 2*ones(size(values.dHPC.decorrelated,1),1) ; ...
             ones(size(values.vHPC.shock,1),1) ; 2*ones(size(values.vHPC.decorrelated,1),1)];

perform2WayUnpairedANOVA(y, structure, condition)
