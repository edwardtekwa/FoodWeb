
function [] = plot_demog_spatial(ptch, Z, B, P, v)
% ptch: patch to plot for time-series
% v   : search rates for all species

scrsz = get(0,'ScreenSize');
%% Zooplankton biomass time series from patch ptch
figure('Position',[1 scrsz(4) scrsz(3)*7/8 scrsz(4)*3/3])
hold on
subplot(3,3,1);
plot(linspace(1,size(Z,3), size(Z,3))-1,log10(squeeze(Z(ptch,1,:))),'linewidth',2,'color',[1 .75 .5]);
xlabel('Time (days)','fontsize',16)
ylabel('Zooplankton biomass (log_{10} g)','fontsize',16)
set(gca,'fontsize',15);
title(['Patch ' num2str(ptch)])
box on;

%% Fish Biomass time series from patch ptch
subplot(3,3,2); hold on;
CM=colormap(jet(128)); % set colormap
caxis([log10(min(P.s.mi)) log10(max(P.s.mi))]); % set colormap range
C = caxis;
for i = 1:P.n
    Ci = fix((log10(P.s.mi(i))-C(1))/(C(2)-C(1))*(size(CM,1)-1))+1;
    plot(linspace(1,size(B,3), size(B,3))-1, log10(squeeze(B(ptch,i,:))),'linewidth',2,'color',CM(Ci,:));
end
xlabel('Time (days)','fontsize',16)
ylabel('Fish biomass (log_{10} g)','fontsize',16)
set(gca,'fontsize',15);
title(['Patch ' num2str(ptch)])
box on;


%% Final biomass size spectra from patch ptch
X = log10(P.s.mi);
Y = log10(B(ptch,:,end));
subplot(3,3,3);
plot(X,Y,'ok','markerfacecolor','k');
set(gca,'fontsize',15);
box on;
xlabel('Fish Body Mass (log_{10} g)','fontsize',16)
ylabel('Fish Biomass (log_{10} g)','fontsize',16)
title(['Patch ' num2str(ptch)])

%% Final distribution of biomasses across all patches
X = 1:P.nx;
Y = log10(B(:,:,end));
subplot(3,3,4:6);
CM=colormap(jet(128)); % set colormap
caxis([log10(min(P.s.mi)) log10(max(P.s.mi))]); % set colormap range
C = caxis;
hold on
for i = 1:P.n
    Ci = fix((log10(P.s.mi(i))-C(1))/(C(2)-C(1))*(size(CM,1)-1))+1;
    plot(X, Y(:,i), 'Color', CM(Ci,:));
end
set(gca,'fontsize',15);
box on;
xlabel('Patch','fontsize',16)
ylabel('Fish Biomass (log_{10} g)','fontsize',16)
at = linspace(log10(min(P.s.mi)), log10(max(P.s.mi)), 5);
b = colorbar('YTick', at);
set(get(b,'ylabel'), 'String', 'Fish Body Mass (log_{10} g)', 'fontsize', 16)

%% Thermal envelopes
X = 1:P.nx;
%Y = squeeze(v); % search rates: spp x patches
Y = bsxfun(@times,squeeze(v),P.s.mi');
subplot(3,3,7:9);
CM=colormap(jet(128)); % set colormap
caxis([log10(min(P.s.mi)) log10(max(P.s.mi))]); % set colormap range
C = caxis;
hold on
for i = 1:P.n
    Ci = fix((log10(P.s.mi(i))-C(1))/(C(2)-C(1))*(size(CM,1)-1))+1;
    plot(X, Y(i,:), 'Color', CM(Ci,:));
end
set(gca,'fontsize',15);
box on;
xlabel('Patch','fontsize',16)
ylabel('Search rate (m^{3} day^{-1})','fontsize',16)


%% Save to file
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 1]);
%print(gcf, '-dpdf', 'Figs/biomass_plots')

end

