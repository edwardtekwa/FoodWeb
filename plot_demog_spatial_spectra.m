
function [] = plot_demog_spatialB(numpatch,Z, B, P, HP,prodZ,prodB,ei)
% ptch: patch to plot for time-series
scrsz = get(0,'ScreenSize');
fig=figure ('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3) scrsz(4)]);

avgWindow=ceil(length(B)*4/5):length(B)-1; %define time window to average long-term statistics

if numpatch==1
    plots=[1];
elseif numpatch==2
    plots=[1 2];
else
    plots=[1 ceil(numpatch/4) ceil(numpatch/2)];
end
%% Fish Biomass time series from patch ptch
for pl=1:length(plots)
    ptch=plots(pl);
%     subplot(3,3,pl); hold on;
%     CM=colormap(jet(128)); % set colormap
%     caxis([log10(min(P.s.mi)) log10(max(P.s.mi))]); % set colormap range
%     C = caxis;
%     plot(linspace(1,size(Z,3), size(Z,3))-1,log10(squeeze(Z(ptch,1,:))),'linewidth',2,'color',[1 .75 .5]);
%     
%     for i = 1:P.n
%         Ci = fix((log10(P.s.mi(i))-C(1))/(C(2)-C(1))*(size(CM,1)-1))+1;
%         plot(linspace(1,size(B,3), size(B,3))-1, log10(squeeze(B(ptch,i,:))),'linewidth',2,'color',CM(Ci,:));
%     end
%     
%     xlabel('Time (days)','fontsize',16)
%     ylabel('Biomass (log_{10} g)','fontsize',16)
%     set(gca,'fontsize',15);
%     title(['Patch ' num2str(ptch)])
%     box on;
%     curYlim=ylim;
    %ylim([log10(eps) curYlim(2)]);
    
    
    %% Average biomass size spectra from patch ptch
    X = log10(P.s.mi);
    %Y = log10(B(ptch,:,end));
    Y = log10(mean(B(ptch,:,avgWindow),3));
    Y2 = log10(mean(prodB(ptch,:,avgWindow),3));
    subplot(1,3,pl);
    %plot(X,Y,'ok','markerfacecolor','b');
    scatter(X,Y,'b');
    hold on
    %plot(X,Y2,'ok','markerfacecolor','r');
    scatter(X,Y2,'r');
%     [hAx,hLine1,hLine2]=plotyy(X,Y,X,Y2);
    hold on
    scatter(log10(P.s.m0),log10(mean(Z(ptch,:,avgWindow),3)),'b');
    scatter(log10(P.s.m0),log10(mean(prodZ(ptch,:,avgWindow),3)),'r');
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
    %ylim([log10(eps) curYlim(2)]);
    
end

%% Average distribution of biomasses across all patches
% X = 1:P.nx;
% %Y = log10(B(:,:,end));
% Y = log10(mean(B(:,:,avgWindow),3));
% subplot(3,3,7:9);
% CM=colormap(jet(128)); % set colormap
% caxis([log10(min(P.s.mi)) log10(max(P.s.mi))]); % set colormap range
% C = caxis;
% hold on
% for i = 1:P.n
%     Ci = fix((log10(P.s.mi(i))-C(1))/(C(2)-C(1))*(size(CM,1)-1))+1;
%     plot(X, Y(:,i), 'Color', CM(Ci,:));
% end
% set(gca,'fontsize',15);
% box on;
% xlabel(['Patch (Basal size=' num2str(P.s.m0) 'g, mean(D)=' num2str(HP.sdm(ei))  ', std(D)=' num2str(HP.sdv(ei)) ', std(optimal T)=' num2str(HP.tv(ei)) ')'],'fontsize',16)
% ylabel('Fish Biomass (log_{10} g)','fontsize',16)
% at = linspace(log10(min(P.s.mi)), log10(max(P.s.mi)), 5);
% b = colorbar('YTick', at);
% set(get(b,'ylabel'), 'String', 'Body Mass (log_{10} g)', 'fontsize', 16)
% curYlim=ylim;
% ylim([log10(eps) curYlim(2)]);


%% Save to file
% set(gcf,'PaperOrientation','landscape');
% set(gcf,'PaperUnits','normalized');
% set(gcf,'PaperPosition', [0 0 1 1]);
% print(gcf, '-dpdf', 'Figs/biomass_plots')

end

