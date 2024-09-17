%%  Aversive
% Members
% Post Aversive
% Post.aversive.members.dHPC(:,sum(isnan(Post.aversive.members.dHPC))>0) = [];
% Post.aversive.members.vHPC(:,sum(isnan(Post.aversive.members.vHPC))>0) = [];

tmp = [];
for i = 1 : size(Post.aversive.members.dHPC,2)
    [ii iii] = min(abs([-0.5 : 0.005 : 0.5]-(-0.05)));
    [ii iiii] = min(abs([-0.5 : 0.005 : 0.5]-0.3));
    
    t = (Post.aversive.members.dHPC(:,i));
    C = (nanmean([t(1:iii) ; t(iiii:end)]));
    C = (nanmean(t(1:iii)));

%     CC = nanstd([t(1:iii) ; t(iiii:end)]);
%     t = (t-C)./CC;
    t = t./C;    
    t = Smooth((t),2,'kernel','gaussian');
    
    tmp = [tmp , t]; clear ii iii iiii C t
end

tmp1 = [];
for i = 1 : size(Pre.aversive.members.dHPC,2)
    [ii iii] = min(abs([-0.5 : 0.005 : 0.5]-(-0.05)));
    [ii iiii] = min(abs([-0.5 : 0.005 : 0.5]-0.3));
    
    t = (Pre.aversive.members.dHPC(:,i));
    C = (nanmean([t(1:iii) ; t(iiii:end)]));
    C = (nanmean(t(1:iii)));
%     CC = nanstd([t(1:iii) ; t(iiii:end)]);
%     t = (t-C)./CC;
    t = t./C;    
    t = Smooth((t),2,'kernel','gaussian');

    tmp1 = [tmp1 , t]; clear ii iii iiii C t
    
end

[ii iii] = min(abs([-0.5 : 0.005 : 0.5]-(-0.2)));
[ii iiii] = min(abs([-0.5 : 0.005 : 0.5]-0.2));

T = nanmean(tmp(iii:iiii,:));
TT = nanmean(tmp1(iii:iiii,:));

t = quantile(nanmean([T ; TT]),[0.25 0.50 0.75]);
% t = quantile(([T - TT]),[0.25 0.50 0.75]);


% T = and([T - TT] < t(3) , [T - TT] > t(2)); %clear t
T = nanmean([T ; TT]) >= t(3);
% T = [T - TT] < t(1);
% T = not([T > TT]);

figure,
subplot(121),plot([-0.5 : 0.005 : 0.5] ,nanmean(tmp(:,T)'),'g'),hold on
ciplot(nanmean(tmp(:,T)')-nansem(tmp(:,T)') , nanmean(tmp(:,T)')+nansem(tmp(:,T)'),[-0.5 : 0.005 : 0.5],'g'),alpha 0.2
plot([-0.5 : 0.005 : 0.5] ,nanmean(tmp1(:,T)'),'b')
ciplot(nanmean(tmp1(:,T)')-nansem(tmp1(:,T)') , nanmean(tmp1(:,T)')+nansem(tmp1(:,T)'),[-0.5 : 0.005 : 0.5],'b'),alpha 0.2
xline(0,'--'),xlim([-0.5 0.5]),ylim([0 6])


PostD = [];
count = 1;
for i = 1 : size(Post.aversive.members.dHPC,2)
    if T(i)
        [ii iii] = min(abs([-0.5 : 0.005 : 0.5]-0.2));
        [ii iiii] = min(abs([-0.5 : 0.005 : 0.5]-0.2));
        
        t = nanmean(tmp(iii:iiii,i));
        
        PostD = [PostD , t]; clear ii iii iiii C t
        count = count+1;
    end
end

PreD = [];
count = 1;
for i = 1 : size(Pre.aversive.members.dHPC,2)
    if T(i)
        [ii iii] = min(abs([-0.5 : 0.005 : 0.5]-0.2));
        [ii iiii] = min(abs([-0.5 : 0.005 : 0.5]-0.2));
        
        t = nanmean(tmp1(iii:iiii,i));
        
        PreD = [PreD , t]; clear ii iii iiii C t
        count = count+1;
    end
end


% Ventral Members
tmp = [];
for i = 1 : size(Post.aversive.members.vHPC,2)
    [ii iii] = min(abs([-0.5 : 0.005 : 0.5]-(-0.05)));
    [ii iiii] = min(abs([-0.5 : 0.005 : 0.5]-0.3));
    
    t = (Post.aversive.members.vHPC(:,i));
    C = (nanmean([t(1:iii) ; t(iiii:end)]));
    C = (nanmean(t(1:iii)));
%     CC = nanstd([t(1:iii) ; t(iiii:end)]);
%     t = (t-C)./CC;
    t = t./C;    
    t = Smooth((t),2,'kernel','gaussian');
    
    tmp = [tmp , t]; clear ii iii iiii C t
end

tmp1 = [];
for i = 1 : size(Pre.aversive.members.vHPC,2)
    [ii iii] = min(abs([-0.5 : 0.005 : 0.5]-(-0.05)));
    [ii iiii] = min(abs([-0.5 : 0.005 : 0.5]-0.3));
    
    t = (Pre.aversive.members.vHPC(:,i));
    C = (nanmean([t(1:iii) ; t(iiii:end)]));
    C = (nanmean(t(1:iii)));
%     CC = nanstd([t(1:iii) ; t(iiii:end)]);
%     t = (t-C)./CC;
    t = t./C;    
    t = Smooth((t),2,'kernel','gaussian');
    
    tmp1 = [tmp1 , t]; clear ii iii iiii C t
    
end

[ii iii] = min(abs([-0.5 : 0.005 : 0.5]-(-0.2)));
[ii iiii] = min(abs([-0.5 : 0.005 : 0.5]-0.2));

T = nanmean(tmp(iii:iiii,:));
TT = nanmean(tmp1(iii:iiii,:));

t = quantile(nanmean([T ; TT]),[0.25 0.5 0.75]);
% t = quantile(([T - TT]),[0.25 0.50 0.75]);

% T = and([T - TT] < t(3) , [T - TT] > t(2)); %clear t
% T = [T - TT] < t(1);
T = nanmean([T ; TT]) >= t(3);
% T = not([T > TT]);

subplot(122),plot([-0.5 : 0.005 : 0.5] ,nanmean(tmp(:,T)'),'g'),hold on
ciplot(nanmean(tmp(:,T)')-nansem(tmp(:,T)') , nanmean(tmp(:,T)')+nansem(tmp(:,T)'),[-0.5 : 0.005 : 0.5],'g'),alpha 0.2
plot([-0.5 : 0.005 : 0.5] ,nanmean(tmp1(:,T)'),'b')
ciplot(nanmean(tmp1(:,T)')-nansem(tmp1(:,T)') , nanmean(tmp1(:,T)')+nansem(tmp1(:,T)'),[-0.5 : 0.005 : 0.5],'b'),alpha 0.2
xline(0,'--'),xlim([-0.5 0.5]),ylim([0 6])


PostV = [];
for i = 1 : size(Post.aversive.members.vHPC,2)
    if T(i)
        [ii iii] = min(abs([-0.5 : 0.005 : 0.5]-0));
        [ii iiii] = min(abs([-0.5 : 0.005 : 0.5]-0.1));
        
        t = nanmean(tmp(iii:iiii,i));
        
        PostV = [PostV , t]; clear ii iii iiii C t
    end
end

PreV = [];
for i = 1 : size(Pre.aversive.members.vHPC,2)
    if T(i)
        [ii iii] = min(abs([-0.5 : 0.005 : 0.5]-0.0));
        [ii iiii] = min(abs([-0.5 : 0.005 : 0.5]-0.1));
        
        t = nanmean(tmp1(iii:iiii,i));
        
        PreV = [PreV , t]; clear ii iii iiii C t
    end
end


figure
subplot(121)
x = [ones(1,size(PreD,2)) , ones(1,size(PostD,2))*2];
y = [PreD , PostD];
scatter(x,y,"filled",'jitter','on', 'jitterAmount',0.1),hold on
scatter([1 2],[nanmean(PreD) , nanmean(PostD)],'filled'),xlim([0 3])
[h p] = signrank(PreD , PostD)

subplot(122)
x = [ones(1,size(PreV,2)) , ones(1,size(PostV,2))*2];
y = [PreV , PostV];
scatter(x,y,"filled",'jitter','on', 'jitterAmount',0.1),hold on
scatter([1 2],[nanmean(PreV) , nanmean(PostV)],'filled'),xlim([0 3])
[h p] = signrank(PreV , PostV)



%%  Reward
% Dorsal Members
tmp = [];
for i = 1 : size(Post.reward.members.dHPC,2)
    [ii iii] = min(abs([-0.5 : 0.005 : 0.5]-(-0.05)));
    [ii iiii] = min(abs([-0.5 : 0.005 : 0.5]-0.15));
    
    t = (Post.reward.members.dHPC(:,i));
    C = (nanmean([t(1:iii) ; t(iiii:end)]));
%     CC = nanstd([t(1:iii) ; t(iiii:end)]);
%     t = (t-C)./CC;
    t = t./C;    
    t = Smooth((t),2,'kernel','gaussian');
    
    tmp = [tmp , t]; clear ii iii iiii C t
end

tmp1 = [];
for i = 1 : size(Pre.reward.members.dHPC,2)
    [ii iii] = min(abs([-0.5 : 0.005 : 0.5]-(-0.05)));
    [ii iiii] = min(abs([-0.5 : 0.005 : 0.5]-0.1));
    
    t = (Pre.reward.members.dHPC(:,i));
    C = (nanmean([t(1:iii) ; t(iiii:end)]));
%     CC = nanstd([t(1:iii) ; t(iiii:end)]);
%     t = (t-C)./CC;
    t = t./C;    
    t = Smooth((t),2,'kernel','gaussian');
    
    tmp1 = [tmp1 , t]; clear ii iii iiii C t
    
end

[ii iii] = min(abs([-0.5 : 0.005 : 0.5]-(-0.0)));
[ii iiii] = min(abs([-0.5 : 0.005 : 0.5]-0.15));

T = nanmean(tmp(iii:iiii,:));
TT = nanmean(tmp1(iii:iiii,:));

t = quantile(T,0.50);

T = T >= t; %clear t

figure,
subplot(121),plot([-0.5 : 0.005 : 0.5] ,nanmean(tmp(:,T)'),'g'),hold on
ciplot(nanmean(tmp(:,T)')-nansem(tmp(:,T)') , nanmean(tmp(:,T)')+nansem(tmp(:,T)'),[-0.5 : 0.005 : 0.5],'g'),alpha 0.2
plot([-0.5 : 0.005 : 0.5] ,nanmean(tmp1(:,T)'),'b')
ciplot(nanmean(tmp1(:,T)')-nansem(tmp1(:,T)') , nanmean(tmp1(:,T)')+nansem(tmp1(:,T)'),[-0.5 : 0.005 : 0.5],'b'),alpha 0.2
xline(0,'--'),xlim([-0.1 0.3]),ylim([0 10])
yline(2,'--')


PostD = [];
count = 1;
for i = 1 : size(Post.reward.members.dHPC,2)
    if T(i)
        [ii iii] = min(abs([-0.5 : 0.005 : 0.5]-0));
        [ii iiii] = min(abs([-0.5 : 0.005 : 0.5]-0.1));
        
        t = max(tmp(iii:iiii,i));
        
        PostD = [PostD , t]; clear ii iii iiii C t
        count = count+1;
    end
end

PreD = [];
count = 1;
for i = 1 : size(Pre.reward.members.dHPC,2)
    if T(i)
        [ii iii] = min(abs([-0.5 : 0.005 : 0.5]-0.0));
        [ii iiii] = min(abs([-0.5 : 0.005 : 0.5]-0.1));
        
        t = max(tmp1(iii:iiii,i));
        
        PreD = [PreD , t]; clear ii iii iiii C t
        count = count+1;
    end
end


% Ventral Members
tmp = [];
for i = 1 : size(Post.reward.members.vHPC,2)
    [ii iii] = min(abs([-0.5 : 0.005 : 0.5]-(-0.05)));
    [ii iiii] = min(abs([-0.5 : 0.005 : 0.5]-0.15));
    
    t = (Post.reward.members.vHPC(:,i));
    C = (nanmean([t(1:iii) ; t(iiii:end)]));
%     CC = nanstd([t(1:iii) ; t(iiii:end)]);
%     t = (t-C)./CC;
    t = t./C;    
    t = Smooth((t),2,'kernel','gaussian');
    
    tmp = [tmp , t]; clear ii iii iiii C t
end

tmp1 = [];
for i = 1 : size(Pre.reward.members.vHPC,2)
    [ii iii] = min(abs([-0.5 : 0.005 : 0.5]-(-0.05)));
    [ii iiii] = min(abs([-0.5 : 0.005 : 0.5]-0.15));
    
    t = (Pre.reward.members.vHPC(:,i));
    C = (nanmean([t(1:iii) ; t(iiii:end)]));
%     CC = nanstd([t(1:iii) ; t(iiii:end)]);
%     t = (t-C)./CC;
    t = t./C;    
    t = Smooth((t),2,'kernel','gaussian');
    
    tmp1 = [tmp1 , t]; clear ii iii iiii C t
    
end

[ii iii] = min(abs([-0.5 : 0.005 : 0.5]-(-0.0)));
[ii iiii] = min(abs([-0.5 : 0.005 : 0.5]-0.1));

T = nanmean(tmp(iii:iiii,:));
TT = nanmean(tmp1(iii:iiii,:));

t = quantile(T,0.50);

T = T >= t; %clear t

subplot(122),plot([-0.5 : 0.005 : 0.5] ,nanmean(tmp(:,T)'),'g'),hold on
ciplot(nanmean(tmp(:,T)')-nansem(tmp(:,T)') , nanmean(tmp(:,T)')+nansem(tmp(:,T)'),[-0.5 : 0.005 : 0.5],'g'),alpha 0.2
plot([-0.5 : 0.005 : 0.5] ,nanmean(tmp1(:,T)'),'b')
ciplot(nanmean(tmp1(:,T)')-nansem(tmp1(:,T)') , nanmean(tmp1(:,T)')+nansem(tmp1(:,T)'),[-0.5 : 0.005 : 0.5],'b'),alpha 0.2
xline(0,'--'),xlim([-0.1 0.3]),ylim([0 10])
yline(2,'--')


PostV = [];
for i = 1 : size(Post.reward.members.vHPC,2)
    if T(i)
        [ii iii] = min(abs([-0.5 : 0.005 : 0.5]-0));
        [ii iiii] = min(abs([-0.5 : 0.005 : 0.5]-0.1));
        
        t = max(tmp(iii:iiii,i));
        
        PostV = [PostV , t]; clear ii iii iiii C t
    end
end

PreV = [];
for i = 1 : size(Pre.reward.members.vHPC,2)
    if T(i)
        [ii iii] = min(abs([-0.5 : 0.005 : 0.5]-0.0));
        [ii iiii] = min(abs([-0.5 : 0.005 : 0.5]-0.1));
        
        t = max(tmp1(iii:iiii,i));
        
        PreV = [PreV , t]; clear ii iii iiii C t
    end
end


figure
subplot(121)
x = [ones(1,size(PreD,2)) , ones(1,size(PostD,2))*2];
y = [PreD , PostD];
scatter(x,y,"filled",'jitter','on', 'jitterAmount',0.1),hold on
scatter([1 2],[nanmean(PreD) , nanmean(PostD)],'filled'),xlim([0 3])
[h p] = signrank(PreD , PostD)

subplot(122)
x = [ones(1,size(PreV,2)) , ones(1,size(PostV,2))*2];
y = [PreV , PostV];
scatter(x,y,"filled",'jitter','on', 'jitterAmount',0.1),hold on
scatter([1 2],[nanmean(PreV) , nanmean(PostV)],'filled'),xlim([0 3])
[h p] = signrank(PreV , PostV)




%% ----- Peak concentration -----
%  Aversive members
tmp = [];
for i = 1 : size(Post.aversive.members.dHPC,2)

    t = (Post.aversive.members.dHPC(:,i));
    t = Smooth((t),2,'kernel','gaussian');
    t = t- min(t);
    t = t./max(t);
    
    tmp = [tmp , t]; clear ii iii iiii C t
end

figure,
subplot(121),plot([-0.5 : 0.005 : 0.5] ,nanmean(tmp'),'g'),hold on
ciplot(nanmean(tmp')-nansem(tmp') , nanmean(tmp')+nansem(tmp'),[-0.5 : 0.005 : 0.5],'g'),alpha 0.2
xline(0,'--'),xlim([-0.1 0.3]),ylim([0 1])

% Ventral Members
tmp = [];
for i = 1 : size(Post.aversive.members.vHPC,2)
    [ii iii] = min(abs([-0.5 : 0.005 : 0.5]-(-0.05)));
    [ii iiii] = min(abs([-0.5 : 0.005 : 0.5]-0.15));
    
    t = (Post.aversive.members.vHPC(:,i));
    t = Smooth((t),2,'kernel','gaussian');
    t = t- min(t);
    t = t./max(t);
    
    tmp = [tmp , t]; clear ii iii iiii C t
end

subplot(122),plot([-0.5 : 0.005 : 0.5] ,nanmean(tmp'),'b'),hold on
ciplot(nanmean(tmp')-nansem(tmp') , nanmean(tmp')+nansem(tmp'),[-0.5 : 0.005 : 0.5],'b'),alpha 0.2
xline(0,'--'),xlim([-0.1 0.3]),ylim([0 1])


%  Reward
% Dorsal members
tmp = [];
for i = 1 : size(Post.reward.members.dHPC,2)
    t = (Post.reward.members.dHPC(:,i));
    t = Smooth((t),2,'kernel','gaussian');
    t = t- min(t);
    t = t./max(t);
    
    tmp = [tmp , t]; clear ii iii iiii C t
end

figure,
subplot(121),plot([-0.5 : 0.005 : 0.5] ,nanmean(tmp'),'g'),hold on
ciplot(nanmean(tmp')-nansem(tmp') , nanmean(tmp')+nansem(tmp'),[-0.5 : 0.005 : 0.5],'g'),alpha 0.2
xline(0,'--'),xlim([-0.1 0.3]),ylim([0 1])


% Ventral Members
tmp = [];
for i = 1 : size(Post.reward.members.vHPC,2)
    
    t = (Post.reward.members.vHPC(:,i));
    t = Smooth((t),2,'kernel','gaussian');
    t = t- min(t);
    t = t./max(t);
    
    tmp = [tmp , t]; clear ii iii iiii C t
end


subplot(122),plot([-0.5 : 0.005 : 0.5] ,nanmean(tmp(:,T)'),'g'),hold on
ciplot(nanmean(tmp(:,T)')-nansem(tmp(:,T)') , nanmean(tmp(:,T)')+nansem(tmp(:,T)'),[-0.5 : 0.005 : 0.5],'g'),alpha 0.2
xline(0,'--'),xlim([-0.1 0.3]),ylim([0 10])
