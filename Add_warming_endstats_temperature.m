%Edward Wong Oct 22, 17
%collect data from Make_warming_MultPtsstats_Parallel0 simulations

clear
set(0,'defaultaxeslinewidth',2)
set(0,'DefaultAxesFontSize',14)

%choose movement scenario to display: (position of movement =[-Inf,-1,0,1,10])
recordYrs=20; %use these number of years at the end of time series to compute mean statistics
%pick files with only one movement rate treatment (meanD)
Cases={'basalSize0.01_meanD3_stdD0'};%'basalSize0.01_meanD-1_stdD0';'basalSize0.01_meanD0_stdD0'}; _pInedible0.5
%Cases={'basalSize0.01_meanD-Inf_stdD0';'basalSize0.01_meanD0_stdD0';'basalSize0.01_meanD1_stdD0';'basalSize0.01_meanD10_stdD0'};
%Cases={'basalSize0.01_meanD-Inf_stdD1';'basalSize0.01_meanD0_stdD1';'basalSize0.01_meanD1_stdD1';'basalSize0.01_meanD10_stdD1'};
%Cases={'basalSize0.01_meanD-Inf_stdD0'};
numCases=length(Cases);
numTemps=2; %include no temp change scenario [0 2 4 6]
zeroLoc=1; %location of zero temperature change when inserted into array of temperature
maxNumIt=80; %maximum number of iterations
%parameters: (record the case parameters at indices in each variable that follows)
meanDispersals=[];
stdDispersals=[];
basalSizes=[];
paramIndices=[];
warmingIndices=[];
%variables for no warming/warming (col1, col2 unless otherwise commented) cases:
%global (for heterotrophs)
totBiomass=[]; %total biomass
totProd=[]; %total production
totRich=[]; %species richness
totBeta=[]; %beta diversity
totTrophicLevel=[]; %mean global trophic level
maxTrophicLevel=[]; %maximum global non-extinct trophic species
totBodyMass=[]; %mean body mass (by biomass)
maxBodyMass=[]; %maximum body mass
maxBiomass=[]; %biomass of most common species
maxProd=[]; %production of most productive species
maxBiomassBS=[]; %body size of most common species
maxProdBS=[]; %body size of most productive species
K_T_ratios=[];
r_T_ratios=[];
transBiomass=[]; %total biomass in transcient period
transProd=[]; %total production in transcient period
consResRatio=[]; %consumer-to-basal resource ratio
fractSpeciesProd=[]; %fraction of species that are net productive
%by space
varBiomass=[];
varProd=[];
varTrophicLevel=[];
%by body mass
allBody_Biomass=[]; %record all non-extinct species. col1: body mass, col2: biomass for no warming cases, col3: biomass for warming cases.
allBody_Prod=[]; %record all non-extinct species. col1: body mass, col2: production for no warming cases, col3: production for warming cases.
BiomBins=[];
ProdBins=[];
BiomLV1Bins=[];
ProdLV1Bins=[];

BS_all=[]; %record parameter values for each replicate
B_all=[]; %record biomasses for each replicate (averaged over last recordYrs)
Prod_all=[]; %record productions for each replicate (averaged over last recordYrs)
Z_all=[]; %record basal biomass for each replicate (averaged over last recordYrs)

originLoc=[]; %location with highest biomass of each species at the end of transcient period
finalLoc=[]; %location with highest biomass of each species at the end
finalLocLV1=[]; %location with highest biomass of each species at the end
rangeShift=[]; %range shift in warming case
Shift_all=[]; %record shift of each species for each replicate
BiomShift=[]; %community biomass shift in warming case


'Select files containing simulation results in first folder'

Files={};
Paths={};

reply1 = 'a';

while (reply1 == 'a')
    [FileName,PathName]=uigetfile('.mat','Simulation data','MultiSelect','on')
    newFiles=cellstr(FileName);
    newPaths=cellstr(PathName);
    Files=[Files, newFiles];
    for f=1:length(newFiles)
        Paths=[Paths, newPaths];
    end
    reply1 = input('Finished? (enter to continue, a to add files, q to quit): ','s');
    if isempty(reply1)
        reply1 = 'continue';
    end
end

numFiles=length(Files);
CaseNumber=1;
%changes
if (reply1 ~='q') %continue
    if CaseNumber==1 %only record chosen movement scenario
        for TempID=1:numTemps %record temperature scenarios
            if TempID<zeroLoc
                TempScenario=TempID;
            elseif TempID==zeroLoc
                TempScenario=0;
            else
                TempScenario=TempID-1;
            end
            iteration=0;
            Positions=contains(Files,Cases{CaseNumber}); %find positions of .mat files that belong to landscape type CaseNumber
            for filePos=1:numFiles
                if(Positions(filePos)==1 && iteration<maxNumIt)
                    iteration=iteration+1;
                    load([Paths{filePos} Files{filePos}]);
                    
                    %add data to aggregate matrices
                    %parameters: (record the case parameters at indices in each variable that follows)
                    findPos1=findstr(Cases{CaseNumber},'basalSize')+length('basalSize');
                    findPos2=findstr(Cases{CaseNumber},'_meanD')-1;
                    findPos3=findstr(Cases{CaseNumber},'_meanD')+length('_meanD');
                    findPos4=findstr(Cases{CaseNumber},'_pInedible')-1;
                    %findPos4=findstr(Cases{CaseNumber},'_stdD')-1;
                    %findPos5=findstr(Cases{CaseNumber},'_stdD')+length('_stdD');
                    basalSizes=[basalSizes;str2num(Cases{CaseNumber}(findPos1:findPos2))];
                    meanDispersals=[meanDispersals;str2num(Cases{CaseNumber}(findPos3:findPos4))];
                    %stdDispersals=[stdDispersals;str2num(Cases{CaseNumber}(findPos5:end))];
                    paramIndices=[paramIndices;TempScenario];
                    warmingIndices=[warmingIndices;[0 1]];
                    
                    %reduce matrices to contain only last point and TempScenario warming
                    BtransEnd=nanmean(Btrans(:,:,end-recordYrs+1:end),3); %biomasses at the end of transcient period
                    BtransEnd(BtransEnd<=eps)=nan;
                    if TempScenario==0
                        B=nanmean(B_yrs(:,:,end-recordYrs+1:end),3); %B_yrs(:,:,end);
                        B(B<eps)=0;
                        gainB=nanmean(gainB_yrs(:,:,end-recordYrs+1:end),3); %gainB_yrs(:,:,end);
                        Z=nanmean(Z_yrs(:,:,end-recordYrs+1:end),3);
                        Z(Z<eps)=0;
                        gainZ=nanmean(gainZ_yrs(:,:,end-recordYrs+1:end),3);
                        TLall=nanmean(TLall_yrs(end-recordYrs+1:end)); %TLall_yrs(end);
                        TLi=nanmean(TLi_yrs(:,:,end-recordYrs+1:end),3);
                        %---------------Change assignments here for different single-species projections---------
                        BLV1=nanmean(BLV4_yrs(:,:,end-recordYrs+1:end),3);
                        gainBLV1=nanmean(gainBLV4_yrs(:,:,end-recordYrs+1:end),3);
                        %----------------------------------------------------------------------------------------
                        Bnan=B;
                        Bnan(Bnan<=eps)=nan;
                        BLV1nan=BLV1;
                        BLV1nan(BLV1nan<=eps)=nan;
                    else
                        B=nanmean(Bw_yrs(:,:,end-recordYrs+1:end,TempScenario),3);
                        B(B<eps)=0;
                        gainB=nanmean(gainBw_yrs(:,:,end-recordYrs+1:end,TempScenario),3);
                        Z=nanmean(Zw_yrs(:,:,end-recordYrs+1:end,TempScenario),3);
                        Z(Z<eps)=0;
                        gainZ=nanmean(gainZw_yrs(:,:,end-recordYrs+1:end,TempScenario),3);
                        TLall=nanmean(TLallw_yrs(end-recordYrs+1:end,TempScenario));
                        TLi=nanmean(TLiw_yrs(:,:,end-recordYrs+1:end,TempScenario),3);
                        %---------------Change assignments here for different single-species projections---------
                        BLV1=nanmean(BLV4w_yrs(:,:,end-recordYrs+1:end,TempScenario),3);
                        gainBLV1=nanmean(gainBLV4w_yrs(:,:,end-recordYrs+1:end,TempScenario),3);
                        %----------------------------------------------------------------------------------------
                        Bnan=B;
                        Bnan(Bnan<=eps)=nan;
                        BLV1nan=BLV1;
                        BLV1nan(BLV1nan<=eps)=nan;
                    end
                    numPatches=size(B,1);
                    
                    
                    %variables for temperature scenarios in full and in single-species projections (col1, col2 unless otherwise commented) cases:
                    %global (for heterotrophs)
                    %totBiomass=[totBiomass; [nansum(B(:)) nansum(BLV1(:))]]; %total biomass
                    %totProd=[totProd; [nansum(gainB(:)) nansum(gainBLV1(:))]]; %total production
                    totBiomass=[totBiomass; [mean(nansum(B,2)) mean(nansum(BLV1,2))]]; %total biomass (as mean biomass/m^3)
                    totProd=[totProd; [mean(nansum(gainB,2)) mean(nansum(gainBLV1,2))]]; %total production (as mean production/(d m^3))
                    totRich=[totRich; [nansum(nansum(B)>eps*numPatches) nansum(nansum(BLV1)>eps*numPatches)]]; %species richness
                    totBeta=[totBeta; [nansum(nansum(B)>eps*numPatches)/nanmean(nansum(B>eps*numPatches,2)) nansum(nansum(BLV1)>eps*numPatches)/nanmean(nansum(BLV1>eps*numPatches,2))]]; %beta diversity
                    totTrophicLevel=[totTrophicLevel; TLall]; %mean global trophic level
                    maxTrophicLevel=[maxTrophicLevel; max(TLi)]; %maximum non-extinct global trophic species
                    totBodyMass=[totBodyMass; [nansum(nansum(B).*P.S(2:end))/nansum(B(:)) nansum(nansum(BLV1).*P.S(2:end))/nansum(BLV1(:)) ]]; %mean body mass (by biomass)
                    tempMaxBodyMass=P.S(max(find(nansum(B)>eps*numPatches))+1);
                    if isempty(tempMaxBodyMass)
                        tempMaxBodyMass=NaN;
                    end
                    tempSingleSpecMaxBodyMass=P.S(max(find(nansum(BLV1)>eps*numPatches))+1);
                    if isempty(tempSingleSpecMaxBodyMass)
                        tempSingleSpecMaxBodyMass=NaN;
                    end
                    maxBodyMass=[maxBodyMass; [tempMaxBodyMass tempSingleSpecMaxBodyMass]]; %maximum body mass
                    [maxB findMaxB]=max(nansum(B));
                    [maxBLV1 findMaxBLV1]=max(nansum(BLV1));
                    [maxP findMaxP]=max(nansum(gainB));
                    [maxPLV1 findMaxPLV1]=max(nansum(gainBLV1));
                    maxBiomass=[maxBiomass; [maxB./nansum(B(:)) maxBLV1./nansum(BLV1(:))]];
                    %maxProd=[maxProd; [maxP./nansum(gainB(:)) maxPLV1./nansum(gainBLV1(:))]];
                    maxProd=[maxProd; [maxP maxPLV1]];
                    maxBiomassBS=[maxBiomassBS; [P.S(findMaxB+1) P.S(findMaxBLV1+1)]];
                    maxProdBS=[maxProdBS; [P.S(findMaxP+1) P.S(findMaxPLV1+1)]];
                    r_T_ratios=[r_T_ratios; r_T_ratio1];
                    K_T_ratios=[K_T_ratios; K_T_ratio1];
                    transBiomass=[transBiomass; sum(sum(mean(Btrans,3)))]; %total biomass in transcient period
                    transProd=[transProd; sum(sum(mean(gainBtrans,3)))]; %total production in transcient period
                    consResRatio=[consResRatio; nansum(B(:))/nansum(Z(:))];
                    fractSpeciesProd=[fractSpeciesProd; [nansum(nansum(gainB)>0)/nansum(nansum(B)>eps*numPatches) nansum(nansum(gainBLV1)>0)/nansum(nansum(BLV1)>eps*numPatches)]];
                
                    %record biomass and production by body-size bins
                    [Ybins,Y2bins] = get_BPbins(B, P, gainB,1);
                    [YLV1bins,Y2LV1bins] = get_BPbins(BLV1, P, gainBLV1,1);
%                     BiomBins=[BiomBins; 10.^Ybins(:,2:end)];
%                     ProdBins=[ProdBins; 10.^Y2bins(:,2:end)];
%                     BiomLV1Bins=[BiomLV1Bins; 10.^YLV1bins(:,2:end)];
%                     ProdLV1Bins=[ProdLV1Bins; 10.^Y2LV1bins(:,2:end)];
                    BiomBins=[BiomBins; Ybins(:,2:end)];
                    ProdBins=[ProdBins; Y2bins(:,2:end)];
                    BiomLV1Bins=[BiomLV1Bins; YLV1bins(:,2:end)];
                    ProdLV1Bins=[ProdLV1Bins; Y2LV1bins(:,2:end)];
                    
                    BS_all=[BS_all; P.s.mi]; %record parameter values for each replicate
                    B_all=cat(3,B_all,B); %record biomasses for each replicate (averaged over last recordYrs)
                    Prod_all=cat(3,Prod_all,gainB); %record productions for each replicate (averaged over last recordYrs)
                    Z_all=[Z_all Z]; %record basal biomass for each replicate
                    
                    %range shift
                    [dummy BtransEndLoc]=max(BtransEnd);
                    BtransEndLoc(isnan(dummy))=nan;
                    originLoc=[originLoc; BtransEndLoc]; %location with highest biomass of each species at the end of transcient period
                    [dummy BLoc]=max(Bnan);
                    BLoc(isnan(dummy))=nan;
                    finalLoc=[finalLoc; BLoc]; %location with highest biomass of each species at the end of warming period
                    [dummy BLV1Loc]=max(BLV1nan);
                    BLV1Loc(isnan(dummy))=nan;
                    finalLocLV1=[finalLocLV1; BLV1Loc]; %location with highest biomass of each species at the end of warming period
                    Shift_all=cat(3,Shift_all,[(BLoc-BtransEndLoc);(BLV1Loc-BtransEndLoc)]);
                    rangeShift=[rangeShift; nanmean(BLoc-BtransEndLoc) nanmean(BLV1Loc-BtransEndLoc)]; %average location shift
                    [dummy originBiomLoc]=max(nansum(BtransEnd,2));
                    [dummy finalBiomLoc]=max(nansum(Bnan,2));
                    [dummy finalBiomLV1Loc]=max(nansum(BLV1nan,2));
                    BiomShift=[BiomShift; finalBiomLoc-originBiomLoc finalBiomLV1Loc-originBiomLoc];
                end
            end
        end
    end
end

numIt=length(paramIndices)/numTemps; %number of iterations per case
tempChanges=P.dT*73000; %temperature changes corresponding to temperature scenario
tempChanges=[tempChanges(1:zeroLoc-1) 0  tempChanges(zeroLoc:end)];

cd = [uint8(cool(5)*255) uint8(ones(5,1))].'; %blue to red gradient color
cd(3,:)=cd(2,:);
cd(2,:)=0;

scrsz = get(0,'ScreenSize');
figs(1)=figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/2 scrsz(4)]);
subplot(4,2,1)
hold on
bl1=boundedline([1:numTemps], nanmean(reshape((totBiomass(:,1)),numIt,[])),[nanstd(reshape((totBiomass(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((totBiomass(:,1)),numIt,[])))).^0.5)]','b','alpha');
boundedline([1:numTemps], nanmean(reshape((totBiomass(:,2)),numIt,[])),[nanstd(reshape((totBiomass(:,2)),numIt,[]))*1.96./((sum(~isnan(reshape((totBiomass(:,2)),numIt,[])))).^0.5)]','--k','alpha','transparency', 0.1);
% p1=plot(sort(reshape(log10(totBiomass(:,1)),numIt,[]))','b'); %plot individual simulation outcomes
% drawnow;
% for i=1:1:length(p1)
%     set(p1(i).Edge, 'ColorBinding','interpolated', 'ColorData',cd)
% end
%set(bl1,'linestyle','none');
drawnow;
set(bl1,'linewidth',2);
set(bl1.Edge, 'ColorBinding','interpolated', 'ColorData',cd);
%scatter([1:numTemps],reshape(log10(totBiomass(:,1)),numIt,[])','b'); %plot individual simulation outcomes
moveTxtPos1=findstr(Cases{1},'meanD')+5;
%moveTxtPos2=findstr(Cases{1},'_pInedible')-1;
%moveTxtPos2=findstr(Cases{1},'_stdD')-1;
title(['connectivity=' Cases{1}(moveTxtPos1) ', n=' num2str(numIt)])
xlabel 'temperature change'
ylabel 'biomass [gm^{-3}]'
xticklabels(tempChanges);
legend off
subplot(4,2,2)
hold on
%gscatter([paramIndices;paramIndices],totProd(:),warmingIndices(:),'br','o')
hold on
bl1=boundedline([1:numTemps], nanmean(reshape((totProd(:,1)),numIt,[])),[nanstd(reshape((totProd(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((totProd(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2); set(bl1.Edge, 'ColorBinding','interpolated', 'ColorData',cd);
boundedline([1:numTemps], nanmean(reshape((totProd(:,2)),numIt,[])),[nanstd(reshape((totProd(:,2)),numIt,[]))*1.96./((sum(~isnan(reshape((totProd(:,2)),numIt,[])))).^0.5)]','--k','alpha','transparency', 0.1);
% p1=plot(sort(reshape((totProd(:,1)),numIt,[]))','b'); %plot individual simulation outcomes
% drawnow;
% for i=1:1:length(p1)
%     set(p1(i).Edge, 'ColorBinding','interpolated', 'ColorData',cd)
% end
%set(bl1,'linestyle','none');
xlabel 'temperature change'
ylabel 'production [gm^{-3}/d]'
xticklabels(tempChanges);
legend off
subplot(4,2,3)
%gscatter([paramIndices;paramIndices],totRich(:),warmingIndices(:),'br','o')
hold on
bl1=boundedline([1:numTemps], nanmean(reshape((totRich(:,1)),numIt,[])),[nanstd(reshape((totRich(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((totRich(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2); set(bl1.Edge, 'ColorBinding','interpolated', 'ColorData',cd);
bl1=boundedline([1:numTemps], nanmean(reshape((totRich(:,2)),numIt,[])),[nanstd(reshape((totRich(:,2)),numIt,[]))*1.96./((sum(~isnan(reshape((totRich(:,2)),numIt,[])))).^0.5)]','--k','alpha','transparency', 0.1);
xlabel 'temperature change'
ylabel 'richness'
xticklabels(tempChanges);
title ''
legend off
subplot(4,2,4)
%gscatter([paramIndices;paramIndices],totBeta(:),warmingIndices(:),'br','o')
hold on
bl1=boundedline([1:numTemps], nanmean(reshape((totBeta(:,1)),numIt,[])),[nanstd(reshape((totBeta(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((totBeta(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2); set(bl1.Edge, 'ColorBinding','interpolated', 'ColorData',cd); 
bl1=boundedline([1:numTemps], nanmean(reshape((totBeta(:,2)),numIt,[])),[nanstd(reshape((totBeta(:,2)),numIt,[]))*1.96./((sum(~isnan(reshape((totBeta(:,2)),numIt,[])))).^0.5)]','--k','alpha','transparency', 0.1);
xlabel 'temperature change'
ylabel 'beta diversity'
xticklabels(tempChanges);
title ''
legend off
subplot(4,2,5)
hold on
bl1=boundedline([1:numTemps], nanmean(reshape((totBodyMass(:,1)),numIt,[])),[nanstd(reshape((totBodyMass(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((totBodyMass(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2); set(bl1.Edge, 'ColorBinding','interpolated', 'ColorData',cd); 
bl1=boundedline([1:numTemps], nanmean(reshape((totBodyMass(:,2)),numIt,[])),[nanstd(reshape((totBodyMass(:,2)),numIt,[]))*1.96./((sum(~isnan(reshape((totBodyMass(:,2)),numIt,[])))).^0.5)]','--k','alpha','transparency', 0.1);
xlabel 'temperature change'
ylabel 'mean log_{10} body size'
xticklabels(tempChanges);
title ''
legend off
subplot(4,2,6)
hold on
bl1=boundedline([1:numTemps], nanmean(reshape((consResRatio(:,1)),numIt,[])),[nanstd(reshape((consResRatio(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((consResRatio(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2); set(bl1.Edge, 'ColorBinding','interpolated', 'ColorData',cd); 
xlabel 'temperature change'
ylabel 'consumer:resource ratio'
% bl1=boundedline([1:numTemps], nanmean(reshape((maxBodyMass(:,1)),numIt,[])),[nanstd(reshape((maxBodyMass(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((maxBodyMass(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2); set(bl1.Edge, 'ColorBinding','interpolated', 'ColorData',cd); 
% bl1=boundedline([1:numTemps], nanmean(reshape((maxBodyMass(:,2)),numIt,[])),[nanstd(reshape((maxBodyMass(:,2)),numIt,[]))*1.96./((sum(~isnan(reshape((maxBodyMass(:,2)),numIt,[])))).^0.5)]','--k','alpha','transparency', 0.1);
% xlabel 'temperature change'
% ylabel 'max log_{10} body size'
xticklabels(tempChanges);
title ''
legend off
subplot(4,2,7)
%gscatter([paramIndices;paramIndices],totTrophicLevel(:),warmingIndices(:),'br','o')
hold on
bl1=boundedline([1:numTemps], nanmean(reshape((totTrophicLevel(:,1)),numIt,[])),[nanstd(reshape((totTrophicLevel(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((totTrophicLevel(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2); set(bl1.Edge, 'ColorBinding','interpolated', 'ColorData',cd); 
xlabel 'temperature change'
ylabel 'mean trophic level'
xticklabels(tempChanges);
title ''
legend off
subplot(4,2,8)
%gscatter([paramIndices;paramIndices],totTrophicLevel(:),warmingIndices(:),'br','o')
hold on
bl1=boundedline([1:numTemps], nanmean(reshape((maxTrophicLevel(:,1)),numIt,[])),[nanstd(reshape((maxTrophicLevel(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((totTrophicLevel(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2); set(bl1.Edge, 'ColorBinding','interpolated', 'ColorData',cd); 
xlabel 'temperature change'
ylabel 'max trophic level'
xticklabels(tempChanges);
title ''
legend off


scrsz = get(0,'ScreenSize');
figs(2)=figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/2 scrsz(4)]);
%max([nanmean(BiomBins(numIt*(zeroLoc-1)+1:numIt*zeroLoc,:)) nanmean(BiomLV1Bins(numIt*(zeroLoc-1)+1:numIt*zeroLoc,:)) nanmean(BiomBins(numIt*(numTemps-1)+1:numIt*numTemps,:)) nanmean(BiomLV1Bins(numIt*(numTemps-1)+1:numIt*numTemps,:))]); 
subplot(4,2,1)
bl1=boundedline([1:size(BiomBins,2)], nanmean(BiomBins(numIt*(zeroLoc-1)+1:numIt*zeroLoc,:)),[nanstd(BiomBins(numIt*(zeroLoc-1)+1:numIt*zeroLoc,:))*1.96./(sum(~isnan(BiomBins(numIt*(zeroLoc-1)+1:numIt*zeroLoc,:))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2); 
hold on
bl1=boundedline([1:size(BiomLV1Bins,2)], nanmean(BiomLV1Bins(numIt*(zeroLoc-1)+1:numIt*zeroLoc,:)),[nanstd(BiomLV1Bins(numIt*(zeroLoc-1)+1:numIt*zeroLoc,:))*1.96./(sum(~isnan(BiomLV1Bins(numIt*(zeroLoc-1)+1:numIt*zeroLoc,:))).^0.5)]','k','alpha');
xticks([1:4:size(BiomBins,2)+1]-0.5);
xticklabels([floor(log10(P.s.mi)):1:ceil(max(log10(P.s.mi)))]);
xlim([4,size(BiomBins,2)+1])
%ylim([-0.05*max([BiomBins(:);BiomLV1Bins(:)]) 0.5*max([BiomBins(:);BiomLV1Bins(:)])])
ylim([0 inf])
%ylim1=ylim;
xlabel 'log_{10}(body size)'
ylabel 'biomass [gm^{-3}]'
title(['connectivity=' Cases{1}(moveTxtPos1) ', n=' num2str(numIt)])

subplot(4,2,2)
bl1=boundedline([1:size(ProdBins,2)], nanmean(ProdBins(numIt*(zeroLoc-1)+1:numIt*zeroLoc,:)),[nanstd(ProdBins(numIt*(zeroLoc-1)+1:numIt*zeroLoc,:))*1.96./(sum(~isnan(ProdBins(numIt*(zeroLoc-1)+1:numIt*zeroLoc,:))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2); 
hold on
bl1=boundedline([1:size(ProdLV1Bins,2)], nanmean(ProdLV1Bins(numIt*(zeroLoc-1)+1:numIt*zeroLoc,:)),[nanstd(ProdLV1Bins(numIt*(zeroLoc-1)+1:numIt*zeroLoc,:))*1.96./(sum(~isnan(ProdLV1Bins(numIt*(zeroLoc-1)+1:numIt*zeroLoc,:))).^0.5)]','k','alpha');
xticks([1:4:size(BiomBins,2)+1]-0.5);
xticklabels([floor(log10(P.s.mi)):1:ceil(max(log10(P.s.mi)))]);
xlim([4,size(BiomBins,2)+1])
%ylim([-0.05*max([ProdBins(:);ProdLV1Bins(:)]) 0.5*max([ProdBins(:);ProdLV1Bins(:)])])
ylim([0 inf])
%ylim2=ylim;
xlabel 'log_{10}(body size)'
ylabel 'production [gm^{-3}/d]'
title 'no warming'

subplot(4,2,3)
bl1=boundedline([1:size(BiomBins,2)], nanmean(BiomBins(numIt*(numTemps-1)+1:numIt*numTemps,:)),[nanstd(BiomBins(numIt*(numTemps-1)+1:numIt*numTemps,:))*1.96./(sum(~isnan(BiomBins(numIt*(numTemps-1)+1:numIt*numTemps,:))).^0.5)]','r','alpha');  drawnow; set(bl1,'linewidth',2);
hold on
bl1=boundedline([1:size(BiomLV1Bins,2)], nanmean(BiomLV1Bins(numIt*(numTemps-1)+1:numIt*numTemps,:)),[nanstd(BiomLV1Bins(numIt*(numTemps-1)+1:numIt*numTemps,:))*1.96./(sum(~isnan(BiomLV1Bins(numIt*(numTemps-1)+1:numIt*numTemps,:))).^0.5)]','k','alpha');
xticks([1:4:size(BiomBins,2)+1]-0.5);
xticklabels([floor(log10(P.s.mi)):1:ceil(max(log10(P.s.mi)))]);
xlim([4,size(BiomBins,2)+1])
%ylim([-0.05*max([BiomBins(:);BiomLV1Bins(:)]) 0.5*max([BiomBins(:);BiomLV1Bins(:)])])
ylim([0 inf])
%ylim3=ylim;
% ylim([0 max([ylim1(2) ylim3(2)])])
ylabel 'biomass [gm^{-3}]'
% subplot(4,2,1)
% ylim([0 max([ylim1(2) ylim3(2)])])

subplot(4,2,4)
bl1=boundedline([1:size(ProdBins,2)], nanmean(ProdBins(numIt*(numTemps-1)+1:numIt*numTemps,:)),[nanstd(ProdBins(numIt*(numTemps-1)+1:numIt*numTemps,:))*1.96./(sum(~isnan(ProdBins(numIt*(numTemps-1)+1:numIt*numTemps,:))).^0.5)]','r','alpha');  drawnow; set(bl1,'linewidth',2);
hold on
bl1=boundedline([1:size(ProdLV1Bins,2)], nanmean(ProdLV1Bins(numIt*(numTemps-1)+1:numIt*numTemps,:)),[nanstd(ProdLV1Bins(numIt*(numTemps-1)+1:numIt*numTemps,:))*1.96./(sum(~isnan(ProdLV1Bins(numIt*(numTemps-1)+1:numIt*numTemps,:))).^0.5)]','k','alpha');
xticks([1:4:size(BiomBins,2)+1]-0.5);
xticklabels([floor(log10(P.s.mi)):1:ceil(max(log10(P.s.mi)))]);
xlim([4,size(BiomBins,2)+1])
%ylim([-0.05*max([ProdBins(:);ProdLV1Bins(:)]) 0.5*max([ProdBins(:);ProdLV1Bins(:)])])
ylim([0 inf])
% ylim4=ylim;
% ylim([0 max([ylim2(2) ylim4(2)])])
xlabel 'log_{10}(body size)'
xlabel 'log_{10}(body size)'
ylabel 'production [gm^{-3}/d]'
title(['warming '  num2str(tempChanges(end)) '^\circC'])
% subplot(4,2,3)
% ylim([0 max([ylim2(2) ylim4(2)])])


subplot(4,2,5)
hold on
bl1=boundedline([1:numTemps], nanmean(reshape((maxBiomass(:,1)),numIt,[])),[nanstd(reshape((maxBiomass(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((maxBiomass(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2); set(bl1.Edge, 'ColorBinding','interpolated', 'ColorData',cd); 
bl1=boundedline([1:numTemps], nanmean(reshape((maxBiomass(:,2)),numIt,[])),[nanstd(reshape((maxBiomass(:,2)),numIt,[]))*1.96./((sum(~isnan(reshape((maxBiomass(:,2)),numIt,[])))).^0.5)]','--k','alpha','transparency', 0.1);
xlabel 'temperature change'
ylabel({'biomass fraction';'from top species'})
xticklabels(tempChanges);
title ''
legend off
subplot(4,2,6)
hold on
%bl1=boundedline([1:numTemps], nanmean(reshape((maxProd(:,1)),numIt,[])),[nanstd(reshape((maxProd(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((maxProd(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2); set(bl1.Edge, 'ColorBinding','interpolated', 'ColorData',cd); 
%bl1=boundedline([1:numTemps], nanmean(reshape((maxProd(:,2)),numIt,[])),[nanstd(reshape((maxProd(:,2)),numIt,[]))*1.96./((sum(~isnan(reshape((maxProd(:,2)),numIt,[])))).^0.5)]','--k','alpha','transparency', 0.1);
bl1=boundedline([1:numTemps], nanmean(reshape((fractSpeciesProd(:,1)),numIt,[])),[nanstd(reshape((fractSpeciesProd(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((fractSpeciesProd(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2); set(bl1.Edge, 'ColorBinding','interpolated', 'ColorData',cd); 
bl1=boundedline([1:numTemps], nanmean(reshape((fractSpeciesProd(:,2)),numIt,[])),[nanstd(reshape((fractSpeciesProd(:,2)),numIt,[]))*1.96./((sum(~isnan(reshape((fractSpeciesProd(:,2)),numIt,[])))).^0.5)]','--k','alpha','transparency', 0.1);
xlabel 'temperature change'
%ylabel 'production from top species'
ylabel({'productive fraction';'of species'})
xticklabels(tempChanges);
ylim([0 1.1])
title ''
legend off
subplot(4,2,7)
hold on
bl1=boundedline([1:numTemps], nanmean(reshape((maxBiomassBS(:,1)),numIt,[])),[nanstd(reshape((maxBiomassBS(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((maxBiomassBS(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2); set(bl1.Edge, 'ColorBinding','interpolated', 'ColorData',cd); 
bl1=boundedline([1:numTemps], nanmean(reshape((maxBiomassBS(:,2)),numIt,[])),[nanstd(reshape((maxBiomassBS(:,2)),numIt,[]))*1.96./((sum(~isnan(reshape((maxBiomassBS(:,2)),numIt,[])))).^0.5)]','--k','alpha','transparency', 0.1);
xlabel 'temperature change'
ylabel({'most common';'body size'})
xticklabels(tempChanges);
title ''
legend off
subplot(4,2,8)
hold on
bl1=boundedline([1:numTemps], nanmean(reshape((maxProdBS(:,1)),numIt,[])),[nanstd(reshape((maxProdBS(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((maxProdBS(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2); set(bl1.Edge, 'ColorBinding','interpolated', 'ColorData',cd); 
bl1=boundedline([1:numTemps], nanmean(reshape((maxProdBS(:,2)),numIt,[])),[nanstd(reshape((maxProdBS(:,2)),numIt,[]))*1.96./((sum(~isnan(reshape((maxProdBS(:,2)),numIt,[])))).^0.5)]','--k','alpha','transparency', 0.1);
xlabel 'temperature change'
ylabel({'most productive';'body size'})
xticklabels(tempChanges);
title ''
legend off

figs(6)=figure ('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/4 scrsz(4)/4]);
bl1=boundedline([1:numTemps], nanmean(reshape(rangeShift(:,1),numIt,[])),[nanstd(reshape(rangeShift(:,1),numIt,[]))*1.96./((sum(~isnan(reshape(rangeShift(:,1),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2); set(bl1.Edge, 'ColorBinding','interpolated', 'ColorData',cd); 
bl1=boundedline([1:numTemps], nanmean(reshape(rangeShift(:,2),numIt,[])),[nanstd(reshape(rangeShift(:,2),numIt,[]))*1.96./((sum(~isnan(reshape(rangeShift(:,2),numIt,[])))).^0.5)]','--k','alpha','transparency', 0.1);
xlabel 'temperature change'
ylabel 'range shift'
title(['connectivity=' Cases{1}(moveTxtPos1) ', n=' num2str(numIt)])
xticklabels(tempChanges);
legend off

plot_demog_spatial_Bavg(B_all,BS_all,Prod_all,Z_all,Shift_all,P,tempChanges,numIt,Cases{1}(moveTxtPos1)); %plot average species contributions over space

% mean(K_T_ratios)
% mean(totBiomass(:,2)./totBiomass(:,1))
% mean(transBiomass./totBiomass(:,1))
% mean(r_T_ratios)
% mean(totProd(:,2)./totProd(:,1))
% mean(transProd./totProd(:,1))