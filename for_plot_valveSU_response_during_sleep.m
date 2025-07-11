%% dHPS SUs
c = curve.id.dHPC;
i = I.dHPC;
tmp = [];
response.shock.dHPC = [];
for ii = 1:size(i,1)
    [~ , iii] = min(abs([-2 : 0.1 : 2]-0));
    [~ , iiii] = min(abs([-2 : 0.1 : 2]-1));
    
    t = i(ii,:);
    tt = and(c(:,1) == t(2) , and(c(:,4) == t(3) , c(:,3) == t(4)));
    tmp = [tmp ; c(tt,:)]; 
    
    C = curve.dHPC(:,tt);
    
    if not(isempty(C))
        response.shock.dHPC = [response.shock.dHPC ; max(C(iii:iiii))];
    end
    
    clear t tt
end

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

c = curve.id.vHPC;
i = I.vHPC;
tmp = [];
response.shock.vHPC = [];
for ii = 1:size(i,1)
    [~ , iii] = min(abs([-2 : 0.1 : 2]-0));
    [~ , iiii] = min(abs([-2 : 0.1 : 2]-1));
    
    t = i(ii,:);
    tt = and(c(:,1) == t(2) , and(c(:,4) == t(3) , c(:,3) == t(4)));
    tmp = [tmp ; c(tt,:)]; 
    
    C = curve.vHPC(:,tt);
    
    if not(isempty(C))
        response.shock.vHPC = [response.shock.vHPC ; max(C(iii:iiii))];
    end
    
    clear t tt
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

%% Plot valve curve
% dHPC
index1 = I.dHPC(:,1)==1;
index2 = not(index1);
c1 = curve.dHPC(:,index1);
m1 = nanmean(c1');
s1 = nansem(c1');
figure,
subplot(121),
plot([-5:0.1:5] , m1,'r'),hold on
ciplot(m1-s1 , m1+s1 , [-5:0.1:5], 'r'), alpha 0.5

c2 = curve.dHPC(:,index2);
m2 = nanmean(c2');
s2 = nansem(c2');
plot([-5:0.1:5] , m2,'k'),hold on
ciplot(m2-s2 , m2+s2 , [-5:0.1:5], 'k'), alpha 0.5

% vHPC
index1 = I.vHPC(:,1)==1;
index2 = not(index1);
c1 = curve.vHPC(:,index1);
m1 = nanmean(c1');
s1 = nansem(c1');
subplot(122),
plot([-5:0.1:5] , m1,'r'),hold on
ciplot(m1-s1 , m1+s1 , [-5:0.1:5], 'r'), alpha 0.5

c2 = curve.vHPC(:,index2);
m2 = nanmean(c2');
s2 = nansem(c2');
plot([-5:0.1:5] , m2,'k'),hold on
ciplot(m2-s2 , m2+s2 , [-5:0.1:5], 'k'), alpha 0.5