%Add_warming_endstats_dispersal_plots.m
%Ed Tekwa Dec 19, 2021

%run Add_warming_endstats_dispersal_fixedDir.m first if collecting new
%data, and comment out the following line. Else the following loads the
%simualtion data from the paper's main results
load('WarmingMovementStats_c1_004 sumLL.mat');

scrsz = get(0,'ScreenSize');
set(0,'defaulttextinterpreter','tex');
set(0, 'defaultAxesTickLabelInterpreter','tex');
set(0, 'defaultLegendInterpreter','tex');
set(0,'defaultaxeslinewidth',2)
set(0,'DefaultAxesFontSize',16)
CM=colormap(jet(128)); % set colormap
%numXPresent=length(moveRates);
numXPresent=6; %6
plotColor='b';

figs(1)=figure ('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/2 scrsz(4)/4]);
subplot(1,2,1)
hold on
bl1=boundedline([1:numCases], nanmean(reshape((totBiomass(:,1)),numIt,[])),[nanstd(reshape((totBiomass(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((totBiomass(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2);
bl1=boundedline([1:numCases], nanmean(reshape((totBiomass(:,3)),numIt,[])),[nanstd(reshape((totBiomass(:,3)),numIt,[]))*1.96./((sum(~isnan(reshape((totBiomass(:,3)),numIt,[])))).^0.5)]','--k','alpha','transparency', 0.1);
xlim([1 numXPresent]); xticks([1:numXPresent]); xticklabels(movementLabels);
xlabel 'movement rate [log(m^2/day)]'
ylabel 'biomass [gm^{-3}]'
legend off
subplot(1,2,2)
hold on
hold on
bl1=boundedline([1:numCases], nanmean(reshape((totProd(:,1)),numIt,[])),[nanstd(reshape((totProd(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((totProd(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2);
bl1=boundedline([1:numCases], nanmean(reshape((totProd(:,3)),numIt,[])),[nanstd(reshape((totProd(:,3)),numIt,[]))*1.96./((sum(~isnan(reshape((totProd(:,3)),numIt,[])))).^0.5)]','--k','alpha','transparency', 0.1);
xlim([1 numXPresent]); xticks([1:numXPresent]); xticklabels(movementLabels);
xlabel 'movement rate [log(m^2/day)]'
ylabel 'production [gm^-{3}/day]'
legend off

%composite graphs (outcomes and effects)
compfig=figure ('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/1.75 scrsz(4)]);
set(compfig,'defaultAxesColorOrder',[[0 0 1]; [1 0 0]]);
subplot(4,2,1)
yyaxis left
hold on
bl1=boundedline([1:numCases], nanmean(reshape((totBiomass(:,1)),numIt,[])),[nanstd(reshape((totBiomass(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((totBiomass(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2);
ylabel 'biomass [gm^{-3}]'
yyaxis right
hold on
pChangeFoodWeb=reshape((totBiomass(:,2)-totBiomass(:,1)),numIt,[]);
pChangeSingleSp=reshape((totBiomass(:,4)-totBiomass(:,3)),numIt,[]);
bl1=boundedline([1:numCases], nanmean(pChangeFoodWeb),[nanstd(pChangeFoodWeb)*1.96./((sum(~isnan(pChangeFoodWeb))).^0.5)]','--r','alpha'); drawnow; set(bl1,'linewidth',2);
bl1=boundedline([1:numCases], nanmean(pChangeSingleSp),[nanstd(pChangeSingleSp)*1.96./((sum(~isnan(pChangeSingleSp))).^0.5)]','--k','alpha');
refl=refline(0,0);
set(refl,'color','r')
xlim([1 numXPresent]); xticks([1:numXPresent]); xticklabels(movementLabels);
title ''
ylabel '\Delta'
legend off
subplot(4,2,2)
yyaxis left
hold on
bl1=boundedline([1:numCases], nanmean(reshape((totProd(:,1)),numIt,[])),[nanstd(reshape((totProd(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((totProd(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2);
ylabel 'production [gm^-{3}/day]'
yyaxis right
hold on
totProd(abs(totProd)>0.1)=NaN; %******set unrealistically large production to NaN
pChangeFoodWeb=reshape((totProd(:,2)-totProd(:,1)),numIt,[]);
pChangeSingleSp=reshape((totProd(:,4)-totProd(:,3)),numIt,[]);
bl1=boundedline([1:numCases], nanmean(pChangeFoodWeb),[nanstd(pChangeFoodWeb)*1.96./((sum(~isnan(pChangeFoodWeb))).^0.5)]','--r','alpha'); drawnow; set(bl1,'linewidth',2);
bl1=boundedline([1:numCases], nanmean(pChangeSingleSp),[nanstd(pChangeSingleSp)*1.96./((sum(~isnan(pChangeSingleSp))).^0.5)]','--k','alpha');
refl=refline(0,0);
set(refl,'color','r')
xlim([1 numXPresent]); xticks([1:numXPresent]); xticklabels(movementLabels);
ylabel '\Delta'
legend off
subplot(4,2,3)
yyaxis left
hold on
bl1=boundedline([1:numCases], nanmean(reshape((totRich(:,1)),numIt,[])),[nanstd(reshape((totRich(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((totRich(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2);
ylabel 'richness'
yyaxis right
hold on
pChangeFoodWeb=reshape((totRich(:,2)-totRich(:,1)),numIt,[]);
pChangeSingleSp=reshape((totRich(:,4)-totRich(:,3)),numIt,[]);
bl1=boundedline([1:numCases], nanmean(pChangeFoodWeb),[nanstd(pChangeFoodWeb)*1.96./((sum(~isnan(pChangeFoodWeb))).^0.5)]','--r','alpha'); drawnow; set(bl1,'linewidth',2);
bl1=boundedline([1:numCases], nanmean(pChangeSingleSp),[nanstd(pChangeSingleSp)*1.96./((sum(~isnan(pChangeSingleSp))).^0.5)]','--k','alpha');
refl=refline(0,0);
set(refl,'color','r')
xlim([1 numXPresent]); xticks([1:numXPresent]); xticklabels(movementLabels);
ylabel '\Delta'
title ''
legend off
subplot(4,2,4)
yyaxis left
hold on
bl1=boundedline([1:numCases], nanmean(reshape((totAlpha(:,1)),numIt,[])),[nanstd(reshape((totAlpha(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((totAlpha(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2);
ylabel 'alpha diversity'
yyaxis right
pChangeFoodWeb=reshape((totAlpha(:,2)-totAlpha(:,1)),numIt,[]);
pChangeSingleSp=reshape((totAlpha(:,4)-totAlpha(:,3)),numIt,[]);
bl1=boundedline([1:numCases], nanmean(pChangeFoodWeb),[nanstd(pChangeFoodWeb)*1.96./((sum(~isnan(pChangeFoodWeb))).^0.5)]','--r','alpha'); drawnow; set(bl1,'linewidth',2);
bl1=boundedline([1:numCases], nanmean(pChangeSingleSp),[nanstd(pChangeSingleSp)*1.96./((sum(~isnan(pChangeSingleSp))).^0.5)]','--k','alpha');
refl=refline(0,0);
set(refl,'color','r')
xlim([1 numXPresent]); xticks([1:numXPresent]); xticklabels(movementLabels);
ylabel '\Delta'
title ''
legend off
subplot(4,2,7)
yyaxis left
hold on
bl1=boundedline([1:numCases], nanmean(reshape((totTrophicLevel(:,1)),numIt,[])),[nanstd(reshape((totTrophicLevel(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((totTrophicLevel(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2);
ylabel 'mean trophic level'
yyaxis right
hold on
pChangeFoodWeb=reshape((totTrophicLevel(:,2)-totTrophicLevel(:,1)),numIt,[]);
bl1=boundedline([1:numCases], nanmean(pChangeFoodWeb),[nanstd(pChangeFoodWeb)*1.96./((sum(~isnan(pChangeFoodWeb))).^0.5)]','--r','alpha'); drawnow; set(bl1,'linewidth',2);
refl=refline(0,0);
set(refl,'color','r')
xlim([1 numXPresent]); xticks([1:numXPresent]); xticklabels(movementLabels);
xlabel 'movement rate [log(m^2/day)]'
ylabel '\Delta'
title ''
legend off
subplot(4,2,8)
yyaxis left
hold on
bl1=boundedline([1:numCases], nanmean(reshape((maxTrophicLevel(:,1)),numIt,[])),[nanstd(reshape((maxTrophicLevel(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((totTrophicLevel(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2);
ylabel 'max trophic level'
yyaxis right
hold on
pChangeFoodWeb=reshape((maxTrophicLevel(:,2)-maxTrophicLevel(:,1)),numIt,[]);
bl1=boundedline([1:numCases], nanmean(pChangeFoodWeb),[nanstd(pChangeFoodWeb)*1.96./((sum(~isnan(pChangeFoodWeb))).^0.5)]','--r','alpha'); drawnow; set(bl1,'linewidth',2);
refl=refline(0,0);
set(refl,'color','r')
xlim([1 numXPresent]); xticks([1:numXPresent]); xticklabels(movementLabels);
xlabel 'movement rate [log(m^2/day)]'
ylabel '\Delta'
title ''
legend off
subplot(4,2,5)
yyaxis left
hold on
bl1=boundedline([1:numCases], nanmean(reshape((totBodyMass(:,1)),numIt,[])),[nanstd(reshape((totBodyMass(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((totBodyMass(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2);
ylabel 'mean log_{10}(body size [g])'
yyaxis right
hold on
pChangeFoodWeb=reshape((totBodyMass(:,2)-totBodyMass(:,1)),numIt,[]);
pChangeSingleSp=reshape((totBodyMass(:,4)-totBodyMass(:,3)),numIt,[]);
bl1=boundedline([1:numCases], nanmean(pChangeFoodWeb),[nanstd(pChangeFoodWeb)*1.96./((sum(~isnan(pChangeFoodWeb))).^0.5)]','--r','alpha'); drawnow; set(bl1,'linewidth',2);
bl1=boundedline([1:numCases], nanmean(pChangeSingleSp),[nanstd(pChangeSingleSp)*1.96./((sum(~isnan(pChangeSingleSp))).^0.5)]','--k','alpha');
refl=refline(0,0);
set(refl,'color','r')
xlim([1 numXPresent]); xticks([1:numXPresent]); xticklabels(movementLabels);
ylabel '\Delta'
title ''
legend off
subplot(4,2,6)
yyaxis left
hold on
bl1=boundedline([1:numCases], nanmean(reshape((consResRatio(:,1)),numIt,[])),[nanstd(reshape((consResRatio(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((consResRatio(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2);
ylabel 'consumer:resource'
yyaxis right
hold on
pChangeFoodWeb=reshape((consResRatio(:,2)-consResRatio(:,1)),numIt,[]);
bl1=boundedline([1:numCases], nanmean(pChangeFoodWeb),[nanstd(pChangeFoodWeb)*1.96./((sum(~isnan(pChangeFoodWeb))).^0.5)]','--r','alpha'); drawnow; set(bl1,'linewidth',2);
refl=refline(0,0);
set(refl,'color','r')
xlim([1 numXPresent]); xticks([1:numXPresent]); xticklabels(movementLabels);
ylabel '\Delta'
title ''
legend off

Assemblagefig=figure ('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/2 scrsz(4)/1.5]);
set(Assemblagefig,'defaultAxesColorOrder',[[0 0 1]; [1 0 0]]);
subplot(2,2,1)
hold on
blmean=boundedline([1:numCases], nanmean(reshape(novelSpatialAssemblage(:,2)*100,numIt,[])),[nanstd(reshape(novelSpatialAssemblage(:,2)*100,numIt,[]))*1.96./((sum(~isnan(reshape(novelSpatialAssemblage(:,2)*100,numIt,[])))).^0.5)]',plotColor,'alpha'); drawnow; set(blmean,'linewidth',2);
bl2=boundedline([1:numCases], nanmean(reshape(novelSpatialAssemblage2(:,2)*100,numIt,[])),[nanstd(reshape(novelSpatialAssemblage2(:,2)*100,numIt,[]))*1.96./((sum(~isnan(reshape(novelSpatialAssemblage2(:,2)*100,numIt,[])))).^0.5)]','-','Cmap',CM(ceil(2*128/6),:),'alpha'); drawnow; set(bl2,'linewidth',1);
bl6=boundedline([1:numCases], nanmean(reshape(novelSpatialAssemblage6(:,2)*100,numIt,[])),[nanstd(reshape(novelSpatialAssemblage6(:,2)*100,numIt,[]))*1.96./((sum(~isnan(reshape(novelSpatialAssemblage6(:,2)*100,numIt,[])))).^0.5)]','-','Cmap',CM(end,:),'alpha'); drawnow; set(bl6,'linewidth',1);
scatter(paramIndices+0.5*(rand(size(paramIndices))-0.5), novelSpatialAssemblage2(:,2)*100,8,CM(ceil(2*128/6),:),'filled');
scatter(paramIndices+0.5*(rand(size(paramIndices))-0.5), novelSpatialAssemblage6(:,2)*100,8,CM(end,:),'filled');

blmeanproj=boundedline([1:numCases], nanmean(reshape(novelSpatialAssemblage(:,4)*100,numIt,[])),[nanstd(reshape(novelSpatialAssemblage(:,4)*100,numIt,[]))*1.96./((sum(~isnan(reshape(novelSpatialAssemblage(:,4)*100,numIt,[])))).^0.5)]','--','Cmap',plotColor,'alpha'); drawnow; set(blmeanproj,'linewidth',2);
bl2proj=boundedline([1:numCases], nanmean(reshape(novelSpatialAssemblage2(:,4)*100,numIt,[])),[nanstd(reshape(novelSpatialAssemblage2(:,4)*100,numIt,[]))*1.96./((sum(~isnan(reshape(novelSpatialAssemblage2(:,4)*100,numIt,[])))).^0.5)]','--','Cmap',CM(ceil(2*128/6),:),'alpha'); drawnow; set(bl2proj,'linewidth',1);
bl6proj=boundedline([1:numCases], nanmean(reshape(novelSpatialAssemblage6(:,4)*100,numIt,[])),[nanstd(reshape(novelSpatialAssemblage6(:,4)*100,numIt,[]))*1.96./((sum(~isnan(reshape(novelSpatialAssemblage6(:,4)*100,numIt,[])))).^0.5)]','--','Cmap',CM(end,:),'alpha'); drawnow; set(bl6proj,'linewidth',1);

ylabel '% species locally novel'
xlim([1 numXPresent]); xticks([1:numXPresent]); xticklabels(movementLabels);

subplot(2,2,3)
hold on
blmean=boundedline([1:numCases], nanmean(reshape(lostSpatialAssemblage(:,2)*100,numIt,[])),[nanstd(reshape(lostSpatialAssemblage(:,2)*100,numIt,[]))*1.96./((sum(~isnan(reshape(lostSpatialAssemblage(:,2)*100,numIt,[])))).^0.5)]',plotColor,'alpha'); drawnow; set(blmean,'linewidth',2);
bl2=boundedline([1:numCases], nanmean(reshape(lostSpatialAssemblage2(:,2)*100,numIt,[])),[nanstd(reshape(lostSpatialAssemblage2(:,2)*100,numIt,[]))*1.96./((sum(~isnan(reshape(lostSpatialAssemblage2(:,2)*100,numIt,[])))).^0.5)]','-','Cmap',CM(ceil(2*128/6),:),'alpha'); drawnow; set(bl2,'linewidth',1);
bl6=boundedline([1:numCases], nanmean(reshape(lostSpatialAssemblage6(:,2)*100,numIt,[])),[nanstd(reshape(lostSpatialAssemblage6(:,2)*100,numIt,[]))*1.96./((sum(~isnan(reshape(lostSpatialAssemblage6(:,2)*100,numIt,[])))).^0.5)]','-','Cmap',CM(end,:),'alpha'); drawnow; set(bl6,'linewidth',1);
scatter(paramIndices+0.5*(rand(size(paramIndices))-0.5), lostSpatialAssemblage2(:,2)*100,8,CM(ceil(2*128/6),:),'filled');
scatter(paramIndices+0.5*(rand(size(paramIndices))-0.5), lostSpatialAssemblage6(:,2)*100,8,CM(end,:),'filled');

blmeanproj=boundedline([1:numCases], nanmean(reshape(lostSpatialAssemblage(:,4)*100,numIt,[])),[nanstd(reshape(lostSpatialAssemblage(:,4)*100,numIt,[]))*1.96./((sum(~isnan(reshape(lostSpatialAssemblage(:,4)*100,numIt,[])))).^0.5)]','--','Cmap',plotColor,'alpha'); drawnow; set(blmeanproj,'linewidth',2);
bl2proj=boundedline([1:numCases], nanmean(reshape(lostSpatialAssemblage2(:,4)*100,numIt,[])),[nanstd(reshape(lostSpatialAssemblage2(:,4)*100,numIt,[]))*1.96./((sum(~isnan(reshape(lostSpatialAssemblage2(:,4)*100,numIt,[])))).^0.5)]','--','Cmap',CM(ceil(2*128/6),:),'alpha'); drawnow; set(bl2proj,'linewidth',1);
bl6proj=boundedline([1:numCases], nanmean(reshape(lostSpatialAssemblage6(:,4)*100,numIt,[])),[nanstd(reshape(lostSpatialAssemblage6(:,4)*100,numIt,[]))*1.96./((sum(~isnan(reshape(lostSpatialAssemblage6(:,4)*100,numIt,[])))).^0.5)]','--','Cmap',CM(end,:),'alpha'); drawnow; set(bl6proj,'linewidth',1);

xlabel 'movement rate [log(m^2/day)]'
ylabel '% species locally extirpated'
xlim([1 numXPresent]); xticks([1:numXPresent]); xticklabels(movementLabels);

subplot(2,2,2)
hold on
blmean=boundedline([1:numCases], nanmean(reshape(novelCoexistence(:,2)*100,numIt,[])),[nanstd(reshape(novelCoexistence(:,2)*100,numIt,[]))*1.96./((sum(~isnan(reshape(novelCoexistence(:,2)*100,numIt,[])))).^0.5)]',plotColor,'alpha'); drawnow; set(blmean,'linewidth',2);
bl22=boundedline([1:numCases], nanmean(reshape(novelCoexistence22(:,2)*100,numIt,[])),[nanstd(reshape(novelCoexistence22(:,2)*100,numIt,[]))*1.96./((sum(~isnan(reshape(novelCoexistence22(:,2)*100,numIt,[])))).^0.5)]','-','Cmap',CM(ceil(2*128/6),:),'alpha'); drawnow; set(bl22,'linewidth',1);
bl26=boundedline([1:numCases], nanmean(reshape(novelCoexistence26(:,2)*100,numIt,[])),[nanstd(reshape(novelCoexistence26(:,2)*100,numIt,[]))*1.96./((sum(~isnan(reshape(novelCoexistence26(:,2)*100,numIt,[])))).^0.5)]','m','alpha'); drawnow; set(bl26,'linewidth',1);
bl66=boundedline([1:numCases], nanmean(reshape(novelCoexistence66(:,2)*100,numIt,[])),[nanstd(reshape(novelCoexistence66(:,2)*100,numIt,[]))*1.96./((sum(~isnan(reshape(novelCoexistence66(:,2)*100,numIt,[])))).^0.5)]','-','Cmap',CM(end,:),'alpha'); drawnow; set(bl66,'linewidth',1);
scatter(paramIndices+0.5*(rand(size(paramIndices))-0.5), novelCoexistence22(:,2)*100,8,CM(ceil(2*128/6),:),'filled');
scatter(paramIndices+0.5*(rand(size(paramIndices))-0.5), novelCoexistence26(:,2)*100,8,'m','filled');
scatter(paramIndices+0.5*(rand(size(paramIndices))-0.5), novelCoexistence66(:,2)*100,8,CM(end,:),'filled');

blmeanproj=boundedline([1:numCases], nanmean(reshape(novelCoexistence(:,4)*100,numIt,[])),[nanstd(reshape(novelCoexistence(:,4)*100,numIt,[]))*1.96./((sum(~isnan(reshape(novelCoexistence(:,4)*100,numIt,[])))).^0.5)]','--','Cmap',plotColor,'alpha'); drawnow; set(blmeanproj,'linewidth',2);
bl22proj=boundedline([1:numCases], nanmean(reshape(novelCoexistence22(:,4)*100,numIt,[])),[nanstd(reshape(novelCoexistence22(:,4)*100,numIt,[]))*1.96./((sum(~isnan(reshape(novelCoexistence22(:,4)*100,numIt,[])))).^0.5)]','--','Cmap',CM(ceil(2*128/6),:),'alpha'); drawnow; set(bl22proj,'linewidth',1);
bl26proj=boundedline([1:numCases], nanmean(reshape(novelCoexistence26(:,4)*100,numIt,[])),[nanstd(reshape(novelCoexistence26(:,4)*100,numIt,[]))*1.96./((sum(~isnan(reshape(novelCoexistence26(:,4)*100,numIt,[])))).^0.5)]','m--','alpha'); drawnow; set(bl26proj,'linewidth',1);
bl66proj=boundedline([1:numCases], nanmean(reshape(novelCoexistence66(:,4)*100,numIt,[])),[nanstd(reshape(novelCoexistence66(:,4)*100,numIt,[]))*1.96./((sum(~isnan(reshape(novelCoexistence66(:,4)*100,numIt,[])))).^0.5)]','--','Cmap',CM(end,:),'alpha'); drawnow; set(bl66proj,'linewidth',1);

ylabel '% coexisting pairs novel'
xlim([1 numXPresent]); xticks([1:numXPresent]); xticklabels(movementLabels);

subplot(2,2,4)
hold on
blmean=boundedline([1:numCases], nanmean(reshape(lostCoexistence(:,2)*100,numIt,[])),[nanstd(reshape(lostCoexistence(:,2)*100,numIt,[]))*1.96./((sum(~isnan(reshape(lostCoexistence(:,2)*100,numIt,[])))).^0.5)]',plotColor,'alpha'); drawnow; set(blmean,'linewidth',2);
bl22=boundedline([1:numCases], nanmean(reshape(lostCoexistence22(:,2)*100,numIt,[])),[nanstd(reshape(lostCoexistence22(:,2)*100,numIt,[]))*1.96./((sum(~isnan(reshape(lostCoexistence22(:,2)*100,numIt,[])))).^0.5)]','-','Cmap',CM(ceil(2*128/6),:),'alpha'); drawnow; set(bl22,'linewidth',1);
bl26=boundedline([1:numCases], nanmean(reshape(lostCoexistence26(:,2)*100,numIt,[])),[nanstd(reshape(lostCoexistence26(:,2)*100,numIt,[]))*1.96./((sum(~isnan(reshape(lostCoexistence26(:,2)*100,numIt,[])))).^0.5)]','m','alpha'); drawnow; set(bl26,'linewidth',1);
bl66=boundedline([1:numCases], nanmean(reshape(lostCoexistence66(:,2)*100,numIt,[])),[nanstd(reshape(lostCoexistence66(:,2)*100,numIt,[]))*1.96./((sum(~isnan(reshape(lostCoexistence66(:,2)*100,numIt,[])))).^0.5)]','-','Cmap',CM(end,:),'alpha'); drawnow; set(bl66,'linewidth',1);
scatter(paramIndices+0.5*(rand(size(paramIndices))-0.5), lostCoexistence22(:,2)*100,8,CM(ceil(2*128/6),:),'filled');
scatter(paramIndices+0.5*(rand(size(paramIndices))-0.5), lostCoexistence26(:,2)*100,8,'m','filled');
scatter(paramIndices+0.5*(rand(size(paramIndices))-0.5), lostCoexistence66(:,2)*100,8,CM(end,:),'filled');

blmeanproj=boundedline([1:numCases], nanmean(reshape(lostCoexistence(:,4)*100,numIt,[])),[nanstd(reshape(lostCoexistence(:,4)*100,numIt,[]))*1.96./((sum(~isnan(reshape(lostCoexistence(:,4)*100,numIt,[])))).^0.5)]','--','Cmap',plotColor,'alpha'); drawnow; set(blmeanproj,'linewidth',2);
bl22proj=boundedline([1:numCases], nanmean(reshape(lostCoexistence22(:,4)*100,numIt,[])),[nanstd(reshape(lostCoexistence22(:,4)*100,numIt,[]))*1.96./((sum(~isnan(reshape(lostCoexistence22(:,4)*100,numIt,[])))).^0.5)]','--','Cmap',CM(ceil(2*128/6),:),'alpha'); drawnow; set(bl22proj,'linewidth',1);
bl26proj=boundedline([1:numCases], nanmean(reshape(lostCoexistence26(:,4)*100,numIt,[])),[nanstd(reshape(lostCoexistence26(:,4)*100,numIt,[]))*1.96./((sum(~isnan(reshape(lostCoexistence26(:,4)*100,numIt,[])))).^0.5)]','m--','alpha'); drawnow; set(bl26proj,'linewidth',1);
bl66proj=boundedline([1:numCases], nanmean(reshape(lostCoexistence66(:,4)*100,numIt,[])),[nanstd(reshape(lostCoexistence66(:,4)*100,numIt,[]))*1.96./((sum(~isnan(reshape(lostCoexistence66(:,4)*100,numIt,[])))).^0.5)]','--','Cmap',CM(end,:),'alpha'); drawnow; set(bl66proj,'linewidth',1);

xlabel 'movement rate [log(m^2/day)]'
ylabel '% coexisting pairs lost'
xlim([1 numXPresent]); xticks([1:numXPresent]); xticklabels(movementLabels);

compfig=figure ('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/1.75 3*scrsz(4)/7]);
set(compfig,'defaultAxesColorOrder',[[0 0 1]; [1 0 0]]);
subplot(2,2,1)
yyaxis left
bl1=boundedline([1:numCases], nanmean(reshape((maxBiomass(:,1)),numIt,[])),[nanstd(reshape((maxBiomass(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((maxBiomass(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2);
ylabel ({'biomass fraction'; 'from top species'})
yyaxis right
hold on
pChangeFoodWeb=reshape((maxBiomass(:,2)-maxBiomass(:,1)),numIt,[]);
pChangeSingleSp=reshape((maxBiomass(:,4)-maxBiomass(:,3)),numIt,[]);
bl1=boundedline([1:numCases], nanmean(pChangeFoodWeb),[nanstd(pChangeFoodWeb)*1.96./((sum(~isnan(pChangeFoodWeb))).^0.5)]','--r','alpha'); drawnow; set(bl1,'linewidth',2);
bl1=boundedline([1:numCases], nanmean(pChangeSingleSp),[nanstd(pChangeSingleSp)*1.96./((sum(~isnan(pChangeSingleSp))).^0.5)]','--k','alpha');
refl=refline(0,0);
set(refl,'color','r')
ylabel '\Delta'
xlim([1 numXPresent]); xticks([1:numXPresent]); xticklabels(movementLabels);
title ''
legend off
subplot(2,2,2)
yyaxis left
bl1=boundedline([1:numCases], nanmean(reshape((fractSpeciesProd(:,1)),numIt,[])),[nanstd(reshape((fractSpeciesProd(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((fractSpeciesProd(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2);
ylabel ({'productive fraction'; 'of species'})
yyaxis right
hold on
pChangeFoodWeb=reshape((fractSpeciesProd(:,2)-fractSpeciesProd(:,1)),numIt,[]);
pChangeSingleSp=reshape((fractSpeciesProd(:,4)-fractSpeciesProd(:,3)),numIt,[]);
bl1=boundedline([1:numCases], nanmean(pChangeFoodWeb),[nanstd(pChangeFoodWeb)*1.96./((sum(~isnan(pChangeFoodWeb))).^0.5)]','--r','alpha'); drawnow; set(bl1,'linewidth',2);
bl1=boundedline([1:numCases], nanmean(pChangeSingleSp),[nanstd(pChangeSingleSp)*1.96./((sum(~isnan(pChangeSingleSp))).^0.5)]','--k','alpha');
refl=refline(0,0);
set(refl,'color','r')
ylabel '\Delta'
xlim([1 numXPresent]); xticks([1:numXPresent]); xticklabels(movementLabels);
legend off
subplot(2,2,3)
yyaxis left
bl1=boundedline([1:numCases], nanmean(reshape((maxBiomassBS(:,1)),numIt,[])),[nanstd(reshape((maxBiomassBS(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((maxBiomassBS(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2);
ylabel ({'most common';'log_{10}(body size [g])'})
yyaxis right
hold on
pChangeFoodWeb=reshape((maxBiomassBS(:,2)-maxBiomassBS(:,1)),numIt,[]);
pChangeSingleSp=reshape((maxBiomassBS(:,4)-maxBiomassBS(:,3)),numIt,[]);
bl1=boundedline([1:numCases], nanmean(pChangeFoodWeb),[nanstd(pChangeFoodWeb)*1.96./((sum(~isnan(pChangeFoodWeb))).^0.5)]','--r','alpha'); drawnow; set(bl1,'linewidth',2);
bl1=boundedline([1:numCases], nanmean(pChangeSingleSp),[nanstd(pChangeSingleSp)*1.96./((sum(~isnan(pChangeSingleSp))).^0.5)]','--k','alpha');
refl=refline(0,0);
set(refl,'color','r')
xlabel 'movement rate [log(m^2/day)]'
ylabel '\Delta'
xlim([1 numXPresent]); xticks([1:numXPresent]); xticklabels(movementLabels);
title ''
legend off
subplot(2,2,4)
yyaxis left
bl1=boundedline([1:numCases], nanmean(reshape((maxProdBS(:,1)),numIt,[])),[nanstd(reshape((maxProdBS(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((maxProdBS(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2);
ylabel ({'most productive';'log_{10}(body size [g])'})
yyaxis right
hold on
pChangeFoodWeb=reshape((maxProdBS(:,2)-maxProdBS(:,1)),numIt,[]);
pChangeSingleSp=reshape((maxProdBS(:,4)-maxProdBS(:,3)),numIt,[]);
bl1=boundedline([1:numCases], nanmean(pChangeFoodWeb),[nanstd(pChangeFoodWeb)*1.96./((sum(~isnan(pChangeFoodWeb))).^0.5)]','--r','alpha'); drawnow; set(bl1,'linewidth',2);
bl1=boundedline([1:numCases], nanmean(pChangeSingleSp),[nanstd(pChangeSingleSp)*1.96./((sum(~isnan(pChangeSingleSp))).^0.5)]','--k','alpha');
refl=refline(0,0);
set(refl,'color','k')
xlabel 'movement rate [log(m^2/day)]'
ylabel '\Delta'
xlim([1 numXPresent]); xticks([1:numXPresent]); xticklabels(movementLabels);
title ''
legend off



shiftfig=figure ('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/2 scrsz(4)/1.5]); %3.5
set(shiftfig,'defaultAxesColorOrder',[[0 0 1]; [1 0 0]]);
MedianShift_move=reshape(All_Shifts(2,:,:)-All_Shifts(1,:,:),numIt*size(All_Shifts,2),[]); %reshape individual species shifts with warming by movement rates
MedianShift_moveLV=reshape(All_Shifts(4,:,:)-All_Shifts(3,:,:),numIt*size(All_Shifts,2),[]); %reshape individual species shifts with warming by movement rates

BS_move=reshape(allBody_Biomass(1,:,:),numIt*size(allBody_Biomass,2),[]); %reshape individual species body sizes by movement rates
OptT_move=reshape(allBody_Biomass(2,:,:),numIt*size(allBody_Biomass,2),[]); %reshape individual species optimal temperature by movement rates
subplot(2,2,1) %centroid shift
hold on
CM=colormap(jet(128)); % set colormap

blmean=boundedline([1:numCases], -nanmean(reshape(CentroidShift(:,1)-CentroidShift0(:,1),numIt,[]))*100/TempChange,[nanstd(reshape(CentroidShift(:,1)-CentroidShift0(:,1),numIt,[]))*(100/TempChange)*1.96./((sum(~isnan(reshape(CentroidShift(:,1)-CentroidShift0(:,1),numIt,[])))).^0.5)]',plotColor,'alpha'); drawnow; set(blmean,'linewidth',2);
bl2=boundedline([1:numCases], -nanmean(reshape(Centroid2Shift(:,1)-Centroid2Shift0(:,1),numIt,[]))*100/TempChange,[nanstd(reshape(Centroid2Shift(:,1)-Centroid2Shift0(:,1),numIt,[]))*(100/TempChange)*1.96./((sum(~isnan(reshape(Centroid2Shift(:,1)-Centroid2Shift0(:,1),numIt,[])))).^0.5)]','-','Cmap',CM(ceil(2*128/6),:),'alpha'); drawnow; set(bl2,'linewidth',1);
scatter(paramIndices+0.5*(rand(size(paramIndices))-0.5), -(Centroid2Shift(:,1)-Centroid2Shift0(:,1))*100/TempChange,8,CM(ceil(2*128/6),:),'filled');
bl6=boundedline([1:numCases], -nanmean(reshape(Centroid6Shift(:,1)-Centroid6Shift0(:,1),numIt,[]))*100/TempChange,[nanstd(reshape(Centroid6Shift(:,1)-Centroid6Shift0(:,1),numIt,[]))*(100/TempChange)*1.96./((sum(~isnan(reshape(Centroid6Shift(:,1)-Centroid6Shift0(:,1),numIt,[])))).^0.5)]','-','Cmap',CM(end,:),'alpha'); drawnow; set(bl6,'linewidth',1);
scatter(paramIndices+0.5*(rand(size(paramIndices))-0.5), -(Centroid6Shift(:,1)-Centroid6Shift0(:,1))*100/TempChange,8,CM(end,:),'filled');

blmeanproj=boundedline([1:numCases], -nanmean(reshape(CentroidShift(:,2)-CentroidShift0(:,2),numIt,[]))*100/TempChange,[nanstd(reshape(CentroidShift(:,2)-CentroidShift0(:,2),numIt,[]))*(100/TempChange)*1.96./((sum(~isnan(reshape(CentroidShift(:,2)-CentroidShift0(:,2),numIt,[])))).^0.5)]','--','Cmap',plotColor,'alpha'); drawnow; set(blmeanproj,'linewidth',2);
bl2proj=boundedline([1:numCases], -nanmean(reshape(Centroid2Shift(:,2)-Centroid2Shift0(:,2),numIt,[]))*100/TempChange,[nanstd(reshape(Centroid2Shift(:,2)-Centroid2Shift0(:,2),numIt,[]))*(100/TempChange)*1.96./((sum(~isnan(reshape(Centroid2Shift(:,2)-Centroid2Shift0(:,2),numIt,[])))).^0.5)]','--','Cmap',CM(ceil(2*128/6),:),'alpha'); drawnow; set(bl2proj,'linewidth',1);
bl6proj=boundedline([1:numCases], -nanmean(reshape(Centroid6Shift(:,2)-Centroid6Shift0(:,2),numIt,[]))*100/TempChange,[nanstd(reshape(Centroid6Shift(:,2)-Centroid6Shift0(:,2),numIt,[]))*(100/TempChange)*1.96./((sum(~isnan(reshape(Centroid6Shift(:,2)-Centroid6Shift0(:,2),numIt,[])))).^0.5)]','--','Cmap',CM(end,:),'alpha'); drawnow; set(bl6proj,'linewidth',1);

refl=refline(0,0);
set(refl,'color','k')
ylabel 'centroid shift %'
ylim([0 100])
xlim([1 numXPresent]); xticks([1:numXPresent]); xticklabels(movementLabels);

subplot(2,2,2) %range contraction
hold on
CM=colormap(jet(128)); % set colormap
blmean=boundedline([1:numCases], -nanmean(reshape(RangeExpansionPerc(:,1)-RangeExpansionPerc0(:,1),numIt,[])),[nanstd(reshape(RangeExpansionPerc(:,1)-RangeExpansionPerc0(:,1),numIt,[]))*1.96./((sum(~isnan(reshape(RangeExpansionPerc(:,1)-RangeExpansionPerc0(:,1),numIt,[])))).^0.5)]',plotColor,'alpha'); drawnow; set(blmean,'linewidth',2);
bl2=boundedline([1:numCases], -nanmean(reshape(Range2ExpansionPerc(:,1)-Range2ExpansionPerc0(:,1),numIt,[])),[nanstd(reshape(Range2ExpansionPerc(:,1)-Range2ExpansionPerc0(:,1),numIt,[]))*1.96./((sum(~isnan(reshape(Range2ExpansionPerc(:,1)-Range2ExpansionPerc0(:,1),numIt,[])))).^0.5)]','-','Cmap',CM(ceil(2*128/6),:),'alpha'); drawnow; set(bl2,'linewidth',1);
scatter(paramIndices+0.5*(rand(size(paramIndices))-0.5), -(Range2ExpansionPerc(:,1)-Range2ExpansionPerc0(:,1)),8,CM(ceil(2*128/6),:),'filled');
bl6=boundedline([1:numCases], -nanmean(reshape(Range6ExpansionPerc(:,1)-Range6ExpansionPerc0(:,1),numIt,[])),[nanstd(reshape(Range6ExpansionPerc(:,1)-Range6ExpansionPerc0(:,1),numIt,[]))*1.96./((sum(~isnan(reshape(Range6ExpansionPerc(:,1)-Range6ExpansionPerc0(:,1),numIt,[])))).^0.5)]','-','Cmap',CM(end,:),'alpha'); drawnow; set(bl6,'linewidth',1);
scatter(paramIndices+0.5*(rand(size(paramIndices))-0.5), -(Range6ExpansionPerc(:,1)-Range6ExpansionPerc0(:,1)),8,CM(end,:),'filled');

blmeanproj=boundedline([1:numCases], -nanmean(reshape(RangeExpansionPerc(:,2)-RangeExpansionPerc0(:,2),numIt,[])),[nanstd(reshape(RangeExpansionPerc(:,2)-RangeExpansionPerc0(:,2),numIt,[]))*1.96./((sum(~isnan(reshape(RangeExpansionPerc(:,2)-RangeExpansionPerc0(:,2),numIt,[])))).^0.5)]','--','Cmap',plotColor,'alpha'); drawnow; set(blmeanproj,'linewidth',2);
bl2proj=boundedline([1:numCases], -nanmean(reshape(Range2ExpansionPerc(:,2)-Range2ExpansionPerc0(:,2),numIt,[])),[nanstd(reshape(Range2ExpansionPerc(:,2)-Range2ExpansionPerc0(:,2),numIt,[]))*1.96./((sum(~isnan(reshape(Range2ExpansionPerc(:,2)-Range2ExpansionPerc0(:,2),numIt,[])))).^0.5)]','--','Cmap',CM(ceil(2*128/6),:),'alpha'); drawnow; set(bl2proj,'linewidth',1);
bl6proj=boundedline([1:numCases], -nanmean(reshape(Range6ExpansionPerc(:,2)-Range6ExpansionPerc0(:,2),numIt,[])),[nanstd(reshape(Range6ExpansionPerc(:,2)-Range6ExpansionPerc0(:,2),numIt,[]))*1.96./((sum(~isnan(reshape(Range6ExpansionPerc(:,2)-Range6ExpansionPerc0(:,2),numIt,[])))).^0.5)]','--','Cmap',CM(end,:),'alpha'); drawnow; set(bl6proj,'linewidth',1);

refl=refline(0,0);
set(refl,'color','k')
ylabel 'range contraction %'
xlim([1 numXPresent]); xticks([1:numXPresent]); xticklabels(movementLabels);

subplot(2,2,3) %leading edge shift
hold on
CM=colormap(jet(128)); % set colormap
blmean=boundedline([1:numCases], -nanmean(reshape(LeadingShift(:,1)-LeadingShift0(:,1),numIt,[]))*100/TempChange,[nanstd(reshape(LeadingShift(:,1)-LeadingShift0(:,1),numIt,[]))*(100/TempChange)*1.96./((sum(~isnan(reshape(LeadingShift(:,1)-LeadingShift0(:,1),numIt,[])))).^0.5)]',plotColor,'alpha'); drawnow; set(blmean,'linewidth',2);
bl2=boundedline([1:numCases], -nanmean(reshape(Leading2Shift(:,1)-Leading2Shift0(:,1),numIt,[]))*100/TempChange,[nanstd(reshape(Leading2Shift(:,1)-Leading2Shift0(:,1),numIt,[]))*(100/TempChange)*1.96./((sum(~isnan(reshape(Leading2Shift(:,1)-Leading2Shift0(:,1),numIt,[])))).^0.5)]','-','Cmap',CM(ceil(2*128/6),:),'alpha'); drawnow; set(bl2,'linewidth',1);
scatter(paramIndices+0.5*(rand(size(paramIndices))-0.5), -(Leading2Shift(:,1)-Leading2Shift0(:,1))*100/TempChange,8,CM(ceil(2*128/6),:),'filled');
bl6=boundedline([1:numCases], -nanmean(reshape(Leading6Shift(:,1)-Leading6Shift0(:,1),numIt,[]))*100/TempChange,[nanstd(reshape(Leading6Shift(:,1)-Leading6Shift0(:,1),numIt,[]))*(100/TempChange)*1.96./((sum(~isnan(reshape(Leading6Shift(:,1)-Leading6Shift0(:,1),numIt,[])))).^0.5)]','-','Cmap',CM(end,:),'alpha'); drawnow; set(bl6,'linewidth',1);
scatter(paramIndices+0.5*(rand(size(paramIndices))-0.5), -(Leading6Shift(:,1)-Leading6Shift0(:,1))*100/TempChange,8,CM(end,:),'filled');

blmeanproj=boundedline([1:numCases], -nanmean(reshape(LeadingShift(:,2)-LeadingShift0(:,2),numIt,[]))*100/TempChange,[nanstd(reshape(LeadingShift(:,2)-LeadingShift0(:,2),numIt,[]))*(100/TempChange)*1.96./((sum(~isnan(reshape(LeadingShift(:,2)-LeadingShift0(:,2),numIt,[])))).^0.5)]','--','Cmap',plotColor,'alpha'); drawnow; set(blmeanproj,'linewidth',2);
bl2proj=boundedline([1:numCases], -nanmean(reshape(Leading2Shift(:,2)-Leading2Shift0(:,2),numIt,[]))*100/TempChange,[nanstd(reshape(Leading2Shift(:,2)-Leading2Shift0(:,2),numIt,[]))*(100/TempChange)*1.96./((sum(~isnan(reshape(Leading2Shift(:,2)-Leading2Shift0(:,2),numIt,[])))).^0.5)]','--','Cmap',CM(ceil(2*128/6),:),'alpha'); drawnow; set(bl2proj,'linewidth',1);
bl6proj=boundedline([1:numCases], -nanmean(reshape(Leading6Shift(:,2)-Leading6Shift0(:,2),numIt,[]))*100/TempChange,[nanstd(reshape(Leading6Shift(:,2)-Leading6Shift0(:,2),numIt,[]))*(100/TempChange)*1.96./((sum(~isnan(reshape(Leading6Shift(:,2)-Leading6Shift0(:,2),numIt,[])))).^0.5)]','--','Cmap',CM(end,:),'alpha'); drawnow; set(bl6proj,'linewidth',1);

refl=refline(0,0);
set(refl,'color','k')
ylabel 'leading edge shift %'
ylim([0 100])
xlim([1 numXPresent]); xticks([1:numXPresent]); xticklabels(movementLabels);
xlabel 'movement rate [log(m^2/day)]'

subplot(2,2,4) %trailing edge shift
hold on
CM=colormap(jet(128)); % set colormap
blmean=boundedline([1:numCases], -nanmean(reshape(TrailingShift(:,1)-TrailingShift0(:,1),numIt,[]))*100/TempChange,[nanstd(reshape(TrailingShift(:,1)-TrailingShift0(:,1),numIt,[]))*(100/TempChange)*1.96./((sum(~isnan(reshape(TrailingShift(:,1)-TrailingShift0(:,1),numIt,[])))).^0.5)]',plotColor,'alpha'); drawnow; set(blmean,'linewidth',2);
bl2mean=boundedline([1:numCases], -nanmean(reshape(Trailing2Shift(:,1)-Trailing2Shift0(:,1),numIt,[]))*100/TempChange,[nanstd(reshape(Trailing2Shift(:,1)-Trailing2Shift0(:,1),numIt,[]))*(100/TempChange)*1.96./((sum(~isnan(reshape(Trailing2Shift(:,1)-Trailing2Shift0(:,1),numIt,[])))).^0.5)]','-','Cmap',CM(ceil(2*128/6),:),'alpha'); drawnow; set(bl2,'linewidth',1);
scatter(paramIndices+0.5*(rand(size(paramIndices))-0.5), -(Trailing2Shift(:,1)-Trailing2Shift0(:,1))*100/TempChange,8,CM(ceil(2*128/6),:),'filled');
bl6mean=boundedline([1:numCases], -nanmean(reshape(Trailing6Shift(:,1)-Trailing6Shift0(:,1),numIt,[]))*100/TempChange,[nanstd(reshape(Trailing6Shift(:,1)-Trailing6Shift0(:,1),numIt,[]))*(100/TempChange)*1.96./((sum(~isnan(reshape(Trailing6Shift(:,1)-Trailing6Shift0(:,1),numIt,[])))).^0.5)]','-','Cmap',CM(end,:),'alpha'); drawnow; set(bl6,'linewidth',1);
scatter(paramIndices+0.5*(rand(size(paramIndices))-0.5), -(Trailing6Shift(:,1)-Trailing6Shift0(:,1))*100/TempChange,8,CM(end,:),'filled');

blmeanproj=boundedline([1:numCases], -nanmean(reshape(TrailingShift(:,2)-TrailingShift0(:,2),numIt,[]))*100/TempChange,[nanstd(reshape(TrailingShift(:,2)-TrailingShift0(:,2),numIt,[]))*(100/TempChange)*1.96./((sum(~isnan(reshape(TrailingShift(:,2)-TrailingShift0(:,2),numIt,[])))).^0.5)]','--','Cmap',plotColor,'alpha'); drawnow; set(blmeanproj,'linewidth',2);
bl2meanproj=boundedline([1:numCases], -nanmean(reshape(Trailing2Shift(:,2)-Trailing2Shift0(:,2),numIt,[]))*100/TempChange,[nanstd(reshape(Trailing2Shift(:,2)-Trailing2Shift0(:,2),numIt,[]))*(100/TempChange)*1.96./((sum(~isnan(reshape(Trailing2Shift(:,2)-Trailing2Shift0(:,2),numIt,[])))).^0.5)]','--','Cmap',CM(ceil(2*128/6),:),'alpha'); drawnow; set(bl2proj,'linewidth',1);
bl6meanproj=boundedline([1:numCases], -nanmean(reshape(Trailing6Shift(:,2)-Trailing6Shift0(:,2),numIt,[]))*100/TempChange,[nanstd(reshape(Trailing6Shift(:,2)-Trailing6Shift0(:,2),numIt,[]))*(100/TempChange)*1.96./((sum(~isnan(reshape(Trailing6Shift(:,2)-Trailing6Shift0(:,2),numIt,[])))).^0.5)]','--','Cmap',CM(end,:),'alpha'); drawnow; set(bl6proj,'linewidth',1);

refl=refline(0,0);
set(refl,'color','k')
ylabel 'trailing edge shift %'
ylim([0 100])
xlim([1 numXPresent]); xticks([1:numXPresent]); xticklabels(movementLabels);
xlabel 'movement rate [log(m^2/day)]'