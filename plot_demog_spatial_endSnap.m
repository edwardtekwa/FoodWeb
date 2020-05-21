%plot_demog_spatial_Btrans.m
%Edward Tekwa March 21, 2019
%plot one simulation

TempChangePos=1;
%run as script:
B=cat(3,Btrans,Bw_yrs(:,:,:,TempChangePos)); %3 refers to position of temperature change (6C warming)
prodB=cat(3,gainBtrans,gainBw_yrs(:,:,:,TempChangePos));

%or use as function:
%function [Gbar,Ybins,Y2bins,avgB,avgP,figs] = plot_demog_spatial_Btrans(B, P, prodB)
%and call it as:
%plot_demog_spatial_Btrans(cat(3,Btrans,Bw_yrs(:,:,:,3)), P, cat(3,gainBtrans,gainBw_yrs(:,:,:,3)))

% ptch: patch to plot for time-series
scrsz = get(0,'ScreenSize');
figs(1)=figure ('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3) scrsz(4)]);

Z=zeros(size(B,1),1,size(B,3));
prodZ=Z;
numpatch=size(B,1);
%avgWindow=ceil(length(B)*1/2):length(B)-1; %define time window to average long-term statistics
%avgWindow=1:length(B);
avgWindow=length(B)-9:length(B); %final (end of warming)
avgWindow1=1:9; %initial (just before warming)
avgWindow2=length(Bw_yrs)-9:length(Bw_yrs); %initial (just before warming)
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
Xb=log10(P.s.m0); %body mass of basal species
%construct histogram bins:
numBins=7; %number of bins for heterotrophs (total number of bins is 1+numBins)
binWidth=(ceil(max(X))-floor(min(X)))/numBins;
xbinLabels=[-2, floor(min(X)):binWidth:ceil(max(X))];
Ybins=zeros(1,numBins+1); %histogram bins for biomass
Y2bins=zeros(1,numBins+1); %histogram bins for productivity

BiomassColor=[0.7 0.7 1];
ProdColor=[1 0.7 0.7];
BiomassErrColor='b';
ProdErrColor='r';

for pl=1:length(plots)
    %% Fish Biomass time series from patch ptch
    ptch=plots(pl);
    subplot(3,3,pl); hold on;
    CM=colormap(jet(128)); % set colormap
    caxis([log10(min(P.s.mi)) log10(max(P.s.mi))]); % set colormap range
    C = caxis;
    %plot(linspace(1,size(Z,3), size(Z,3))-1,log10(squeeze(Z(ptch,1,:))),'linewidth',2,'color','k');
    
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
    
    %% Average biomass size spectra from patch ptch
    %Y = log10(B(ptch,:,end));
    Y = log10(mean(B(ptch,:,avgWindow),3));
    Y2 = log10(mean(prodB(ptch,:,avgWindow),3));
    subplot(3,3,3+pl);
    %plot(X,Y,'ok','markerfacecolor','b');
    %plot(X,Y2,'ok','markerfacecolor','r');
    plot(X,Y2,'-','Color',ProdErrColor,'LineWidth',2);
    hold on
    %s=scatter(X,Y2,[],CM(fix((X-C(1))/(C(2)-C(1))*(size(CM,1)-1))+1,:),'filled');
    %set(s,'MarkerEdgeColor',ProdErrColor);
    plot(X,Y,'-','Color',BiomassErrColor,'LineWidth',2);
    s=scatter(X,Y,[],CM(fix((X-C(1))/(C(2)-C(1))*(size(CM,1)-1))+1,:),'filled');
    set(s,'MarkerEdgeColor',BiomassErrColor);
    %     [hAx,hLine1,hLine2]=plotyy(X,Y,X,Y2);
%    s=scatter(log10(P.s.m0),log10(mean(Z(ptch,:,avgWindow),3)),'k','filled');
%    set(s,'MarkerEdgeColor',BiomassErrColor);
%    s=scatter(log10(P.s.m0),log10(mean(prodZ(ptch,:,avgWindow),3)),'k','filled');
%    set(s,'MarkerEdgeColor',ProdErrColor);
    %     hLine1.LineWidth=2;
    %     hLine1.Color='b';
    %     hLine2.LineWidth=2;
    %     hLine2.Color='r';
    set(gca,'fontsize',15);
    box on;
    xlabel('Body Mass (log_{10} g)','fontsize',16)
    ylabel({['\fontsize{16} {\color{blue}Biomass (log_{10} g) '],...
        ['\color{red}Productivity (log_{10} g/day)}']})
    %     ylabel(hAx(1),'Fish Biomass (log_{10} g)','fontsize',16)
    %     ylabel(hAx(2),'Fish Productivity (log_{10} g)','fontsize',16)
    title(['Patch ' num2str(ptch)])
    curYlim=ylim;
    %xlim([floor(log10(P.s.m0)) ceil(max(X))]);
    xlim([floor(min(X)) ceil(max(X))]);
    ylim([-35 1]);
    
    %biomass and productivity histograms
    subplot(3,3,6+pl);
    Yptch=mean(B(ptch,:,avgWindow),3); %all heterotroph species biomasses
    Y2ptch=mean(prodB(ptch,:,avgWindow),3); %all heteroph species productivities
%    Ybptch=log10(mean(Z(ptch,:,avgWindow),3)); %basal species biomass
%    Y2bptch=log10(mean(prodZ(ptch,:,avgWindow),3)); %basal species productivity
    
%    Ybins(1)=Ybptch; %first bin is basal species
%    Y2bins(1)=Y2bptch;
%    sortedTimeBiomass=sort(log10(Z(ptch,1,avgWindow)));
%    sortedTimeProd=sort(log10(prodZ(ptch,1,avgWindow)));
%     YbinsTimeHi(1)=sortedTimeBiomass(:,:,ceil(lengthWindow*(1-0.025))); %top 2.5% of biomass in time
%     YbinsTimeLo(1)=sortedTimeBiomass(:,:,ceil(lengthWindow*(0.025))); %bottom 2.5% of biomass in time
%     Y2binsTimeHi(1)=sortedTimeProd(:,:,ceil(lengthWindow*(1-0.025))); %top 2.5% of biomass in time
%     Y2binsTimeLo(1)=sortedTimeProd(:,:,ceil(lengthWindow*(0.025))); %bottom 2.5% of biomass in time
%     
    %YbinsTimestd(1)=std(log10(Z(ptch,1,avgWindow))); %(vertical) standard deviation of log biomass over time (averaged over space)
    %Y2binsTimeStd(1)=std(log10(prodZ(ptch,1,avgWindow))); %(vertical) standard deviation of log biomass over time (averaged over space)
    
    for bin=2:numBins+1 %other bins are heterotrophs
        lowerX=floor(min(X))+(bin-2)*binWidth;
        upperX=floor(min(X))+(bin-1)*binWidth;
        Ybins(bin)=log10(sum(Yptch(find(X>=lowerX & X<upperX))));
        Y2bins(bin)=log10(sum(Y2ptch(find(X>=lowerX & X<upperX))));
        sortedTimeBiomass=sort(log10(sum(B(ptch,find(X>=lowerX & X<upperX),avgWindow))));
        sortedTimeProd=sort(log10(sum(prodB(ptch,find(X>=lowerX & X<upperX),avgWindow))));
        YbinsTimeHi(bin)=sortedTimeBiomass(:,:,ceil(lengthWindow*(1-0.025))); %top 2.5% of biomass in time
        YbinsTimeLo(bin)=sortedTimeBiomass(:,:,ceil(lengthWindow*(0.025))); %bottom 2.5% of biomass in time
        Y2binsTimeHi(bin)=sortedTimeProd(:,:,ceil(lengthWindow*(1-0.025))); %top 2.5% of biomass in time
        Y2binsTimeLo(bin)=sortedTimeProd(:,:,ceil(lengthWindow*(0.025))); %bottom 2.5% of biomass in time
        %YbinsTimeStd(bin)=std(log10(sum(B(ptch,find(X>=lowerX & X<upperX),avgWindow))));
        %Y2binsTimeStd(bin)=std(log10(sum(prodB(ptch,find(X>=lowerX & X<upperX),avgWindow))));
    end
    
    b=bar(Ybins-floor(log10(eps)),'FaceColor',BiomassColor);
    hold on
    b=bar(Y2bins-floor(log10(eps)),'FaceColor',ProdColor);
    errorbar([1:1+numBins]-0.04,Ybins-floor(log10(eps)),Ybins-YbinsTimeLo,YbinsTimeHi-Ybins,'.','Color',BiomassErrColor,'LineWidth',2);
    errorbar([1:1+numBins]+0.04,Y2bins-floor(log10(eps)),Y2bins-Y2binsTimeLo,Y2binsTimeHi-Y2bins,'.','Color',ProdErrColor,'LineWidth',2);
    %errorbar([1:8]-0.01,Ybins-floor(log10(eps)),YbinsTimeStd,'.k','LineWidth',2);
    %errorbar([1:8]+0.01,Y2bins-floor(log10(eps)),Y2binsTimeStd,'.k','LineWidth',2);
    xlabel('Body Mass (log_{10} g)','fontsize',16)
    ylabel({['\fontsize{16} {\color{blue}Biomass (log_{10} g) '],...
        ['\color{red}Productivity (log_{10} g/day)}']});
    title(['Patch ' num2str(ptch)])
    ax=gca;
    ax.XTick=[1:2+numBins]-0.5;
    ax.XTickLabel=xbinLabels;
    ax.YTick=[0:2:20];
    ax.YTickLabel=[floor(log10(eps)):2:16-floor(log10(eps))];
    xlim([1.5 numBins+1.5]);
    ylim([16+log10(eps) 18]);
end


%% Average distribution of biomasses across all patches
figs(2)=figure ('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3) scrsz(4)]);
Xp = 1:P.nx;
%Y = log10(B(:,:,end));
Y = log10(mean(B(:,:,avgWindow),3));
CM=colormap(jet(128)); % set colormap
caxis([log10(min(P.s.mi)) log10(max(P.s.mi))]); % set colormap range
C = caxis;

subplot(3,3,4:6);
%yyaxis left
%[plot2,line1,line2]=plotyy(Xp,0,Xp,P.T);
hold on
%set(line2,'LineWidth',2,'Color','r')
%ylim(plot2(1),[0 5])
SpeciesDistr=plot(([(mean(Bw_yrs(:,:,avgWindow1),3))]),'LineWidth',2);
%SpeciesDistr(1).FaceColor='k';
for i = 1:P.n
    Ci = fix((log10(P.s.mi(i))-C(1))/(C(2)-C(1))*(size(CM,1)-1))+1;
    SpeciesDistr(i).Color=CM(Ci,:);
    %plot(Xp, Y(:,i), 'Color', CM(Ci,:),'linewidth',2);
end
set(gca,'fontsize',15);
box on;
xlabel('Patch','fontsize',16)
ylabel('Biomass (g)','fontsize',16)
title('Initial species distributions')
at = linspace(floor(log10(min(P.s.mi))), ceil(log10(max(P.s.mi))), 8);
b = colorbar('YTick',[0:1/7:1],'YTickLabel',at);
set(get(b,'ylabel'), 'String', 'Body Mass (log_{10} g)', 'fontsize', 16)
hold on
%yyaxis right
%plot(Xp,P.T,'linewidth',2);
%ylabel('Temperature','fontsize',16)
xlim([0.8 numpatch+0.2]);

subplot(3,3,7:9);
%plot(Xp, log10(mean(Z(:,:,avgWindow),3)), 'Color', 'k','linewidth',2);
SpeciesDistr=plot((mean(Bw_yrs(:,:,avgWindow2),3)),'LineWidth',2);
hold on
for i = 1:P.n
    Ci = fix((log10(P.s.mi(i))-C(1))/(C(2)-C(1))*(size(CM,1)-1))+1;
    SpeciesDistr(i).Color=CM(Ci,:);
    %plot(Xp, Y(:,i), 'Color', CM(Ci,:),'linewidth',2);
end
set(gca,'fontsize',15);
box on;
%xlabel(['Patch (Basal size=' num2str(P.s.m0) 'g, mean(D)=' num2str(HP.sdm(ei))  ', std(D)=' num2str(HP.sdv(ei)) ', std(optimal T)=' num2str(HP.tv(ei)) ')'],'fontsize',16)
xlabel(['Patch (Basal size=' num2str(P.s.m0) 'g'],'fontsize',16)
ylabel('Biomass (g)','fontsize',16)
title('Final species distributions','fontsize',16)
at = linspace(floor(log10(min(P.s.mi))), ceil(log10(max(P.s.mi))), 8);
b = colorbar('YTick',[0:1/7:1],'YTickLabel',at);
set(get(b,'ylabel'), 'String', 'Body Mass (log_{10} g)', 'fontsize', 16)
curYlim=ylim;
%ylim([0 2.5]);
xlim([0.8 numpatch+0.2]);

%yyaxis right
AlphaRichnessDistr=sort(sum((B(:,:,avgWindow))>eps,2),3);
meanAlphaRichness=mean(AlphaRichnessDistr,3);
loAlphaRichness=AlphaRichnessDistr(:,:,ceil(lengthWindow*0.025));
hiAlphaRichness=AlphaRichnessDistr(:,:,ceil(lengthWindow*0.975));
%bar(meanAlphaRichness,'FaceColor','none','LineWidth',2,'EdgeColor',[0 0.45 0.74]);
%hold on
%er=errorbar([1:numpatch],meanAlphaRichness,meanAlphaRichness-loAlphaRichness,hiAlphaRichness-meanAlphaRichness,'.');
%set(er,'LineWidth',2);
%ylabel('Local richness (\alpha)')
% ax=gca;
% ax.YTick=[0.1 1 10];
% ax.YTickLabel=[-1 0 1];
% ylim([0.1 10]);

%histograms (10 bins) of average global biomasses, productivities, and their
%coefficients of variation
%figure ('Color', [1 1 1])
subplot(3,3,1)
% Yall=(mean(mean(B(:,:,avgWindow),3))); %all heterotroph species biomasses
% Y2all=(mean(mean(prodB(:,:,avgWindow),3))); %all heteroph species productivities
Yball=log10(mean(mean(Z(:,:,avgWindow),3))); %basal species biomass
Y2ball=log10(mean(mean(prodZ(:,:,avgWindow),3))); %basal species productivity

Ybins(1)=Yball; %first bin is basal species
Y2bins(1)=Y2ball;
% YbinsTimeStd(1)=(std(mean(Z(:,1,avgWindow)))); %(vertical) log of standard deviation of biomass over time (averaged over space)
% YbinsSpaceStd(1)=(mean(std(Z(:,1,avgWindow)))); %(horizontal) log of standard deviation of biomass over space (averaged over time)
% Y2binsTimeStd(1)=(std(mean(prodZ(:,1,avgWindow)))); %(vertical) log of standard deviation of biomass over time (averaged over space)
% Y2binsSpaceStd(1)=(mean(std(prodZ(:,1,avgWindow)))); %(horizontal) standard deviation of log biomass over space (averaged over time)
sortedTimeBiomass=sort(log10(mean(Z(:,1,avgWindow))));
sortedTimeProd=sort(log10(mean(prodZ(:,1,avgWindow))));
YbinsTimeHi(1)=sortedTimeBiomass(:,:,ceil(lengthWindow*(1-0.025))); %top 2.5% of biomass in time
YbinsTimeLo(1)=sortedTimeBiomass(:,:,ceil(lengthWindow*(0.025))); %bottom 2.5% of biomass in time
Y2binsTimeHi(1)=sortedTimeProd(:,:,ceil(lengthWindow*(1-0.025))); %top 2.5% of biomass in time
Y2binsTimeLo(1)=sortedTimeProd(:,:,ceil(lengthWindow*(0.025))); %bottom 2.5% of biomass in time

sortedSpaceBiomass=sort(log10(mean(Z(:,1,avgWindow),3)));
sortedSpaceProd=sort(log10(mean(prodZ(:,1,avgWindow),3)));
YbinsSpaceHi(1)=sortedSpaceBiomass(ceil(numpatch*(1-0.025))); %top 2.5% of biomass in time
YbinsSpaceLo(1)=sortedSpaceBiomass(ceil(numpatch*(0.025))); %bottom 2.5% of biomass in time
Y2binsSpaceHi(1)=sortedSpaceProd(ceil(numpatch*(1-0.025))); %top 2.5% of biomass in time
Y2binsSpaceLo(1)=sortedSpaceProd(ceil(numpatch*(0.025))); %bottom 2.5% of biomass in time


for bin=2:numBins+1 %other bins are heterotrophs
    lowerX=floor(min(X))+(bin-2)*binWidth;
    upperX=floor(min(X))+(bin-1)*binWidth;
    %     Ybins(bin)=log10(sum(Yall(find(X>=lowerX & X<upperX))));
    %     Y2bins(bin)=log10(sum(Y2all(find(X>=lowerX & X<upperX))));
    sortedTimeBiomass=sort(log10(mean(sum(B(:,find(X>=lowerX & X<upperX),avgWindow),2))));
    sortedTimeProd=sort(log10(mean(sum(prodB(:,find(X>=lowerX & X<upperX),avgWindow),2))));
    Ybins(bin)=mean(sortedTimeBiomass);
    Y2bins(bin)=mean(sortedTimeProd);
    YbinsTimeHi(bin)=sortedTimeBiomass(:,:,ceil(lengthWindow*(1-0.025))); %top 2.5% of biomass in time
    YbinsTimeLo(bin)=sortedTimeBiomass(:,:,ceil(lengthWindow*(0.025))); %bottom 2.5% of biomass in time
    Y2binsTimeHi(bin)=sortedTimeProd(:,:,ceil(lengthWindow*(1-0.025))); %top 2.5% of biomass in time
    Y2binsTimeLo(bin)=sortedTimeProd(:,:,ceil(lengthWindow*(0.025))); %bottom 2.5% of biomass in time
    %     YbinsTimeStd(bin)=log10(std(sum(mean(B(:,find(X>=lowerX & X<upperX),avgWindow)))));
    %     YbinsSpaceStd(bin)=log10(mean(std(sum(B(:,find(X>=lowerX & X<upperX),avgWindow),2))));
    %     Y2binsTimeStd(bin)=log10(std(sum(mean(prodB(:,find(X>=lowerX & X<upperX),avgWindow)))));
    %     Y2binsSpaceStd(bin)=log10(mean(std(sum(prodB(:,find(X>=lowerX & X<upperX),avgWindow),2))));
    sortedSpaceBiomass=sort(log10(mean(sum(B(:,find(X>=lowerX & X<upperX),avgWindow),2),3)));
    sortedSpaceProd=sort(log10(mean(sum(prodB(:,find(X>=lowerX & X<upperX),avgWindow),2),3)));
    YbinsSpaceHi(bin)=sortedSpaceBiomass(ceil(numpatch*(1-0.025))); %top 2.5% of biomass in time
    YbinsSpaceLo(bin)=sortedSpaceBiomass(ceil(numpatch*(0.025))); %bottom 2.5% of biomass in time
    Y2binsSpaceHi(bin)=sortedSpaceProd(ceil(numpatch*(1-0.025))); %top 2.5% of biomass in time
    Y2binsSpaceLo(bin)=sortedSpaceProd(ceil(numpatch*(0.025))); %bottom 2.5% of biomass in time
    
end
horErrScale=0.05;
%subplot(1,2,1)
b=bar(Ybins-floor(log10(eps)),'FaceColor',BiomassColor);
hold on
b=bar(Y2bins-floor(log10(eps)),'FaceColor',ProdColor);
errorbar([1:1+numBins]-0.04,Ybins-floor(log10(eps)),Ybins-YbinsTimeLo,YbinsTimeHi-Ybins,'.','Color',BiomassErrColor,'LineWidth',2);
%errorbar([1:1+numBins]-0.01,Ybins-floor(log10(eps)),YbinsTimeStd,'.k','LineWidth',2);
%her=herrorbar([1:1+numBins],Ybins-floor(log10(eps)),YbinsSpaceStd/10,'.k');
her=herrorbar([1:1+numBins],Ybins-floor(log10(eps)),(Ybins-YbinsSpaceLo)*horErrScale,(YbinsSpaceHi-Ybins)*horErrScale,'.');
set(her,'LineWidth',2,'Color',BiomassErrColor);
errorbar([1:1+numBins]+0.04,Y2bins-floor(log10(eps)),Y2bins-Y2binsTimeLo,Y2binsTimeHi-Y2bins,'.','Color',ProdErrColor,'LineWidth',2);
%errorbar([1:1+numBins]+0.01,Y2bins-floor(log10(eps)),Y2binsTimeStd,'.k','LineWidth',2);
%her2=herrorbar([1:1+numBins],Y2bins-floor(log10(eps)),Y2binsSpaceStd/10,'.k');
her2=herrorbar([1:1+numBins],Y2bins-floor(log10(eps)),(Y2bins-Y2binsSpaceLo)*horErrScale,(Y2binsSpaceHi-Y2bins)*horErrScale,'.');
set(her2,'LineWidth',2,'Color',ProdErrColor);
xlabel('Body Mass (log_{10} g)','fontsize',16)
ylabel({['\fontsize{16} {\color{blue}Biomass (log_{10} g) '],...
    ['\color{red}Productivity (log_{10} g/day)}']});
title(['Global means, movement rate=' num2str(log10(P.diff(1)))])
ax=gca;
ax.XTick=[1:2+numBins]-0.5;
ax.XTickLabel=xbinLabels;
ax.YTick=[0:2:18];
ax.YTickLabel=[floor(log10(eps)):2:16-floor(log10(eps))];
xlim([0.5 numBins+1.5]);
ylim([16+log10(eps) 18]);

%subplot(3,3,2)
%bar plot with columns representing, in sequence: global
%1. heterotroph biomass
%2. biomass temporal var
%3. biomass spatial var
%4. productivity
%5. prod temp var
%6. prod spatial var
%7. gamma richness
%8. beta richness
%9. gamma Shannon's div
%10. beta Shannon's div
%11. most common species (body mass)
%12. most common species (biomass)

GlobalBiomasses=sort(sum(mean(B(:,:,avgWindow),1),2),3); %mean global heterotroph biomass
Gbar(1)=mean(GlobalBiomasses);
Glo(1)=GlobalBiomasses(:,:,ceil(lengthWindow*(0.025)));
Ghi(1)=GlobalBiomasses(:,:,ceil(lengthWindow*(0.975)));

GlobalProds=sort(sum(mean(prodB(:,:,avgWindow),1),2),3); %mean global heterotroph biomass
Gbar(2)=mean(GlobalProds);
Glo(2)=GlobalProds(:,:,ceil(lengthWindow*(0.025)));
Ghi(2)=GlobalProds(:,:,ceil(lengthWindow*(0.975)));

GammaRichnessDistr=sort(sum(sum(B(:,:,avgWindow),1)>eps,2),3);
Gbar(3)=mean(GammaRichnessDistr);
Glo(3)=GammaRichnessDistr(:,:,ceil(lengthWindow*(0.025)));
Ghi(3)=GammaRichnessDistr(:,:,ceil(lengthWindow*(0.975)));

GammaRichnessDistrUnsort=sum(sum(B(:,:,avgWindow),1)>eps,2);
MeanAlphaRichnessDistrUnsort=mean(sum((B(:,:,avgWindow))>eps,2));
BetaRichnessDistr=sort(GammaRichnessDistrUnsort./MeanAlphaRichnessDistrUnsort,3);
Gbar(4)=mean(BetaRichnessDistr);
Glo(4)=BetaRichnessDistr(:,:,ceil(lengthWindow*(0.025)));
Ghi(4)=BetaRichnessDistr(:,:,ceil(lengthWindow*(0.975)));

% subplot(1,2,2)
% b=bar(Y2bins-floor(log10(eps)),'k');
% xlabel('Body Mass (log_{10} g)','fontsize',16)
% ylabel('Productivity (log_{10} g/day)')
% ax=gca;
% ax.XTick=[1:8]-0.5;
% ax.XTickLabel=xbinLabels;
% ax.YTick=[0:2:18];
% ax.YTickLabel=[floor(log10(eps)):2:16-floor(log10(eps))];
% ylim([12 17]);

%plot alpha & gamma diversities and catch rates

%% Save to file
% set(gcf,'PaperOrientation','landscape');
% set(gcf,'PaperUnits','normalized');
% set(gcf,'PaperPosition', [0 0 1 1]);
% print(gcf, '-dpdf', 'Figs/biomass_plots')
avgB=[(mean(Z(:,:,avgWindow),3)),(mean(B(:,:,avgWindow),3))];
avgP=[(mean(prodZ(:,:,avgWindow),3)),(mean(prodB(:,:,avgWindow),3))];

%savefig(figs,[FoodWebFile '_' num2str(i) '_TimeSeries_noWarming'],'compact');
%     close(figs);

