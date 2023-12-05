
% Marchenko Pastur Distribution
% Ref :
% Marchenko,V. A., Pastur, L. A. (1967) "Distribution of eigenvalues for some sets of
% random matrices", Mat. Sb. (N.S.), 72(114):4, 507�536
%
% --- INPUTS ---
%
% SpksTrains: Matrix of Spikes (rows: cells / columns: time bins)
%
% --- OUTPUTS ---
%
% m: float, surrogated eigen value using Marchenko Pastur Distribution
%  Morici Juan Facundo


function m = marchenko(SpksTrains)


N=size(SpksTrains,2); %number of cells
T=size(SpksTrains,1); %time bins
% Ratio of matrix dimensions
c=N/T;
% Sample
%x=rand(N,T);  % Uniform distribution
s=std(SpksTrains(:));
% spectral matrix
r=SpksTrains'*SpksTrains/T;
%eigenvalues
l=eig(r);
% Probability Density Function 
% number of points for measurement.
n=50;
% Boundaries 
a=(s^2)*(1-sqrt(c))^2;
b=(s^2)*(1+sqrt(c))^2;
[f,lambda]=hist(l,linspace(a,b,n));
% Normalization
f=f/sum(f);
% Theoretical pdf
ft=@(lambda,a,b,c) (1./(2*pi*lambda*c*s^(2))).*sqrt((b-lambda).*(lambda-a));
F=ft(lambda,a,b,c);
% Processing numerical pdf
F=F/sum(F);
F(isnan(F))=0;
% % Results
% figure;
% h=bar(lambda,f);
% set(h,'FaceColor',[.75 .75 .8]);
% set(h,'LineWidth',0.25);
% xlabel('Eigenvalue \lambda');
% ylabel(' Probability Density Function f(\lambda)');
% title(' Marchenko鳳astur distribution');
% lmin=min(l);
% lmax=max(l);
% %axis([-1 2*lmax 0 max(f)+max(f)/4]);
% hold on;
% plot(lambda,F,'g','LineWidth',2);
% [~,mm] = max(F); 
% m = lambda(mm);
m = [a,b];
% hold off;
end
