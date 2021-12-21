%plot_demog_spatial_end.m
%Edward Tekwa March 21, 2019
%plot one simulation

TempChangePos=1;
%run as script:
B=cat(3,Btrans,Bw_yrs(:,:,:,TempChangePos));

% ptch: patch to plot for time-series
scrsz = get(0,'ScreenSize');
figs(1)=figure ('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3) scrsz(4)]);

Z=zeros(size(B,1),1,size(B,3));
prodZ=Z;
numpatch=size(B,1);
avgWindow=length(B)-9:length(B); %final (end of warming)
avgWindow1=length(B)-209:length(B)-200; %initial (just before warming)
avgWindow2=length(B)-9:length(B); %final
lengthWindow=length(avgWindow);

if numpatch==1
    plots=[1];
elseif numpatch==2
    plots=[1 2];
else
    plots=[2 ceil(numpatch/2) numpatch-1];
end

%prepare global biomass and productivity histograms
X = log10(P.s.mi);
Xb=log10(P.s.m0); %body size of basal species

BiomassColor=[0.7 0.7 1];
ProdColor=[1 0.7 0.7];
BiomassErrColor='b';
ProdErrColor='r';

for pl=1:length(plots)
    %% Fish Biomass time series from patch ptch
    ptch=plots(pl);
    subplot(3,3,pl); hold on;
    CM=colormap(jet(128)); % set colormap
    caxis([round(log10(min(P.s.mi))) round(log10(max(P.s.mi)))]); % set colormap range
    C = caxis;
    for i = 1:P.n
        Ci = fix((log10(P.s.mi(i))-C(1))/(C(2)-C(1))*(size(CM,1)-1))+1;
        plot(linspace(1,size(B,3), size(B,3))-1, log10(squeeze(B(ptch,i,:))),'linewidth',2,'color',CM(Ci,:));
    end
    xlabel('time (yr)','fontsize',16)
    ylabel('log_{10}(biomass [gm^{-3}])','fontsize',16)
    set(gca,'fontsize',15);
    title(['patch ' num2str(ptch) ', T=' num2str(P.T(ptch))])
    box on;
    curYlim=ylim;
    ylim([log10(eps) 2]);
    line([800 800],[log10(eps) 2],'LineWidth',2,'Color','k')
    
  
end


subplot(3,3,4:6);
hold on
SpeciesDistr=plot(([(mean(B(:,:,avgWindow1),3))]),'LineWidth',2);
for i = 1:P.n
    Ci = fix((log10(P.s.mi(i))-C(1))/(C(2)-C(1))*(size(CM,1)-1))+1;
    SpeciesDistr(i).Color=CM(Ci,:);
end
set(gca,'fontsize',15);
box on;
xlabel('Patch','fontsize',16)
ylabel('Biomass (g)','fontsize',16)
title('Initial species distributions')
at = linspace(floor(log10(min(P.s.mi))), ceil(log10(max(P.s.mi))), 7);
b = colorbar('YTick',[0:1/6:1],'YTickLabel',at);
set(get(b,'ylabel'), 'String', 'body size (log_{10}g)', 'fontsize', 16)
hold on
xlim([0.8 numpatch+0.2]);


subplot(3,3,7:9);
SpeciesDistr=plot((mean(B(:,:,avgWindow2),3)),'LineWidth',2);
hold on
for i = 1:P.n
    Ci = fix((log10(P.s.mi(i))-C(1))/(C(2)-C(1))*(size(CM,1)-1))+1;
    SpeciesDistr(i).Color=CM(Ci,:);
end
set(gca,'fontsize',15);
box on;
xlabel(['Patch (Basal size=' num2str(P.s.m0) 'g'],'fontsize',16)
ylabel('Biomass (g)','fontsize',16)
title('Final species distributions','fontsize',16)
b = colorbar('YTick',[0:1/6:1],'YTickLabel',at);
set(get(b,'ylabel'), 'String', 'body size (log_{10}g)', 'fontsize', 16)
curYlim=ylim;
xlim([0.8 numpatch+0.2]);

AlphaRichnessDistr=sort(sum((B(:,:,avgWindow))>eps,2),3);
meanAlphaRichness=mean(AlphaRichnessDistr,3);
loAlphaRichness=AlphaRichnessDistr(:,:,ceil(lengthWindow*0.025));
hiAlphaRichness=AlphaRichnessDistr(:,:,ceil(lengthWindow*0.975));
