%plotEnsembles
%Edward Tekwa 5/3/17
%plot ensemble data from 3 scenarios
set(0,'DefaultFigureVisible', 'on'); %use this to suppress plot display but still runs through the stats
    
scrsz = get(0,'ScreenSize');

%first, load ensemble simulations with no movement
load ('Food web 12-Apr-2017 18:21:29.mat')
EnsembleGbar1=EnsembleGbar;
EnsembleGBbins1=EnsembleGBbins;
EnsembleprodBbins1=EnsembleprodBbins;
EnsembleGBbins1(EnsembleGBbins1==-Inf)=log10(eps);
EnsembleprodBbins1(EnsembleprodBbins1==-Inf)=log10(eps);
avgB1=avgB;
avgP1=avgP;
Ps1=Ps;
for i=1:length(Ps)
    size1(:,i)=Ps(i).S;
end
fig=figure ('Color', [1 1 1],'Position',[1 scrsz(2)/2 scrsz(3)/2 scrsz(4)/2]);
[S,AX,BigAx,H,HAx]=plotmatrix(EnsembleGbar1);
title ('no movement')
xlabel(AX(4,1),'total biomass')
xlabel(AX(4,2),'total productivity')
xlabel(AX(4,3),'richness')
xlabel(AX(4,4),'beta diversity')
ylabel(AX(1,1),'total biomass')
ylabel(AX(2,1),'total productivity')
ylabel(AX(3,1),'richness')
ylabel(AX(4,1),'beta diversity')


%second, load ensemble simulations with movement but no variation in
%movement
load ('Food web 14-Apr-2017 18:08:02.mat')
EnsembleGbar2=EnsembleGbar;
EnsembleGBbins2=EnsembleGBbins;
EnsembleprodBbins2=EnsembleprodBbins;
EnsembleGBbins2(EnsembleGBbins2==-Inf)=log10(eps);
EnsembleprodBbins2(EnsembleprodBbins2==-Inf)=log10(eps);
avgB2=avgB;
avgP2=avgP;
Ps2=Ps;
for i=1:length(Ps)
    size2(:,i)=Ps(i).S;
end
fig=figure ('Color', [1 1 1],'Position',[1 scrsz(2)/2 scrsz(3)/2 scrsz(4)/2]);
[S,AX,BigAx,H,HAx]=plotmatrix(EnsembleGbar2);
title ('with movement')
xlabel(AX(4,1),'total biomass')
xlabel(AX(4,2),'total productivity')
xlabel(AX(4,3),'richness')
xlabel(AX(4,4),'beta diversity')
ylabel(AX(1,1),'total biomass')
ylabel(AX(2,1),'total productivity')
ylabel(AX(3,1),'richness')
ylabel(AX(4,1),'beta diversity')

%third, load ensemble simulations with movement and variation in movement
load ('Food web 27-Apr-2017 22:18:03.mat')
EnsembleGbar3=EnsembleGbar;
EnsembleGBbins3=EnsembleGBbins;
EnsembleprodBbins3=EnsembleprodBbins;
EnsembleGBbins3(EnsembleGBbins3==-Inf)=log10(eps);
EnsembleprodBbins3(EnsembleprodBbins3==-Inf)=log10(eps);
avgB3=avgB;
avgP3=avgP;
Ps3=Ps;
for i=1:length(Ps)
    size3(:,i)=Ps(i).S;
end
fig=figure ('Color', [1 1 1],'Position',[1 scrsz(2)/2 scrsz(3)/2 scrsz(4)/2]);
[S,AX,BigAx,H,HAx]=plotmatrix(EnsembleGbar3);
title ('with movement and movement variation')
xlabel(AX(4,1),'total biomass')
xlabel(AX(4,2),'total productivity')
xlabel(AX(4,3),'richness')
xlabel(AX(4,4),'beta diversity')
ylabel(AX(1,1),'total biomass')
ylabel(AX(2,1),'total productivity')
ylabel(AX(3,1),'richness')
ylabel(AX(4,1),'beta diversity')

%compute means and variances of body mass weighted by biomass or
%productivity within each simulation instance
meanB1=sum((reshape(sum(avgB1(:,2:101,:),1),100,20)./(10.^size1(2:101,:))))./EnsembleGbar1(:,1)';
meanB2=sum((reshape(sum(avgB2(:,2:101,:),1),100,20)./(10.^size2(2:101,:))))./EnsembleGbar2(:,1)';
meanB3=sum((reshape(sum(avgB3(:,2:101,:),1),100,20)./(10.^size3(2:101,:))))./EnsembleGbar3(:,1)';


fig=figure ('Color', [1 1 1],'Position',[1 scrsz(2)/2 scrsz(3)/2 scrsz(4)/2]);
subplot(2,2,1)
violin([EnsembleGbar1(:,1),EnsembleGbar2(:,1),EnsembleGbar3(:,1)],'medc',[]);
legend off
ylabel 'biomass'
subplot(2,2,2)
violin([EnsembleGbar1(:,2),EnsembleGbar2(:,2),EnsembleGbar3(:,2)],'medc',[]);
legend off
ylabel 'productivity'
subplot(2,2,3)
violin([EnsembleGbar1(:,3),EnsembleGbar2(:,3),EnsembleGbar3(:,3)],'medc',[]);
legend off
ylabel 'richness'
subplot(2,2,4)
violin([EnsembleGbar1(:,4),EnsembleGbar2(:,4),EnsembleGbar3(:,4)],'medc',[]);
legend off
ylabel 'beta diversity'

fig=figure ('Color', [1 1 1],'Position',[1 scrsz(2)/2 scrsz(3)/2 scrsz(4)/2]);
% scatter(reshape(meshgrid(1:8,1:20),[],1),reshape(EnsembleGBbins1,[],1),'or','jitter','on')
% hold on
% scatter(reshape(meshgrid(1:8,1:20),[],1),reshape(EnsembleGBbins2,[],1),'ob','jitter','on')
% scatter(reshape(meshgrid(1:8,1:20),[],1),reshape(EnsembleGBbins3,[],1),'og','jitter','on')
subplot(2,2,1)
bar([mean(10.^(EnsembleGBbins1(:,2:8)),1); mean(10.^(EnsembleGBbins2(:,2:8)),1); mean(10.^(EnsembleGBbins3(:,2:8)),1) ],'stack')
xticks=[1 2 3];
xticklabels={'no movement' 'with movement' 'movement+variation'};
ylabel 'total biomass'
title 'body mass distribution'
subplot(2,2,2)
bar((10.^(EnsembleGBbins1(:,2:8))),'stack')
title 'without movement'
xlabel 'simulation #'
ylabel 'total biomass'
xlim([0 21])
subplot(2,2,3)
bar((10.^(EnsembleGBbins2(:,2:8))),'stack')
title 'with movement'
xlabel 'simulation #'
ylabel 'total biomass'
xlim([0 21])
subplot(2,2,4)
bar((10.^(EnsembleGBbins3(:,2:8))),'stack')
title 'movement+variation'
xlabel 'simulation #'
ylabel 'total biomass'
xlim([0 21])
