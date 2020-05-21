
%Edward Wong Oct 22, 17
%collect data from Make_warming_MultPtsstats_Parallel0 simulations

%clear
set(0,'defaultaxeslinewidth',2)
set(0,'DefaultAxesFontSize',14)

TimeData=string(datetime);

%choose temperature change scenario to display: (position of tempChange =[-2 2 4 6])
TempScenario=3; %position is without 0 change position [+2 +4 +6]
recordYrs=20; %use these number of years at the end of time series to compute mean statistics
%Cases={'basalSize0.01_meanD-Inf_stdD0';'basalSize0.01_meanD0_stdD0';'basalSize0.01_meanD3_stdD0';'basalSize0.01_meanD6_stdD0';'basalSize0.01_meanD9_stdD0';'basalSize0.01_meanD12_stdD0'};
Cases={'basalSize0.01_meanD-Inf';'basalSize0.01_meanD0';'basalSize0.01_meanD3';'basalSize0.01_meanD6';'basalSize0.01_meanD9';'basalSize0.01_meanD12'};
%Cases={'basalSize0.001_meanD-Inf_stdD0';'basalSize0.001_meanD0_stdD0';'basalSize0.001_meanD3_stdD0';'basalSize0.001_meanD6_stdD0';'basalSize0.001_meanD9_stdD0';'basalSize0.001_meanD12_stdD0'};
numCases=length(Cases);
moveRates=[-Inf,0,3,6,9,12];
%moveRates=[-Inf,0,3,6];
caseIts=[]; %record number of iterations for each case

%quantile fits:
BPquantiles=[];
BPpercDiff=[];
freeparstart=[0.7 0.7];
freeparmin=[0.4 0.4];
freeparmax=[0.9 0.9];

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
totAlpha=[]; %average local alpha diversity
totBeta=[]; %beta diversity
totTrophicLevel=[]; %mean global trophic level
maxTrophicLevel=[]; %maximum global non-extinct trophic species
totBodyMass=[]; %mean body mass (by biomass)
maxBodyMass=[]; %maximum body mass
maxBiomass=[]; %biomass of most common species
maxProd=[]; %production of most productive species
maxBiomassBS=[]; %body size of most common species
maxProdBS=[]; %body size of most productive species
consResRatio=[]; %consumer-to-basal resource ratio
fractSpeciesProd=[]; %fraction of species that are net productive
%by space
varBiomass=[];
varProd=[];
varTrophicLevel=[];
%by body mass
allBody_Biomass=[]; %record all non-extinct species. col1: body mass, col2: biomass for no warming cases, col3: biomass for warming cases.
allBody_Prod=[]; %record all non-extinct species. col1: body mass, col2: production for no warming cases, col3: production for warming cases.
%range shift
originLoc=[]; %location with highest biomass of each species at the end of transcient period
finalLoc=[]; %location with highest biomass of each species at the end
finalLocLV1=[]; %location with highest biomass of each species at the end
finalwLoc=[]; %location with highest biomass of each species at the end of warming period
finalwLocLV1=[]; %location with highest biomass of each species at the end of warming period
rangeShift0=[]; %range shift in no warming case
rangeShift=[]; %range shift in warming case
BiomShift0=[]; %community biomass shift in no warming case
BiomShift=[]; %community biomass shift in warming case
InitCentroid=[]; %mean location of each species at the end of transcient period
FinalCentroid=[]; %mean location of each species at the end
FinalLVCentroid=[]; %projected mean location of each species at the end
FinalwCentroid=[]; %mean location of each species at the end of warming period
FinalwLVCentroid=[]; %mean location of each species at the end of warming period
CentroidShift0=[]; %centroid shift with no warming
CentroidShift=[]; %centroid shift with warming
InitTrailing=[];
InitLeading=[];
FinalTrailing=[];
FinalLeading=[];
FinalwTrailing=[];
FinalwLeading=[];
FinalLVTrailing=[];
FinalLVLeading=[];
FinalwLVTrailing=[];
FinalwLVLeading=[];
TrailingShift0=[]; %average trailing edge shift under no warming
TrailingShift=[]; %average trailing edge shift under warming
LeadingShift0=[]; %average trailing edge shift under no warming
LeadingShift=[]; %average trailing edge shift under warming


%species-specific variables
All_BodySize=[]; %record all surviving species body sizes in all scenarios
All_OptTemp=[]; %optimal temperature
All_InitTemp=[]; %initial patch temperature
All_Move=[]; %movement rate
All_TempChange=[]; %temperature change over 200 years
All_ShiftLag=[]; %lag in species shift
All_Shifts=[]; %all median location shifts
All_CentroidShifts=[]; %all mean (centroid) shifts
All_TrailingShifts=[];
All_LeadingShifts=[];


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

if (reply1 ~='q') %continue
    for CaseNumber=1:numCases %first, check for minimum number of iterations for each case
        iteration=0;
        Positions=contains(Files,Cases{CaseNumber}); %find positions of .mat files that belong to movement rate treatment of CaseNumber
        for filePos=1:numFiles
            if(Positions(filePos)==1)
                iteration=iteration+1;
                caseIts(CaseNumber)=iteration; %update number of iterations for each case
            end
        end
    end
    numIt=min(caseIts); %minimum number of iterations for all cases
    
    for CaseNumber=1:numCases
        iteration=0;
        Positions=contains(Files,Cases{CaseNumber}); %find positions of .mat files that belong to landscape type CaseNumber
        for filePos=1:numFiles
            if(Positions(filePos)==1 && iteration<numIt)
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
                paramIndices=[paramIndices;CaseNumber];
                warmingIndices=[warmingIndices;[0 1]];
                
                %reduce matrices to contain only last point and TempScenario warming
                B=nanmean(B_yrs(:,:,end-recordYrs+1:end),3); %B_yrs(:,:,end);
                gainB=nanmean(gainB_yrs(:,:,end-recordYrs+1:end),3); %gainB_yrs(:,:,end);
                Z=nanmean(Z_yrs(:,:,end-recordYrs+1:end),3);
                gainZ=nanmean(gainZ_yrs(:,:,end-recordYrs+1:end),3);
                TLall=nanmean(TLall_yrs(end-recordYrs+1:end)); %TLall_yrs(end);
                TLi=nanmean(TLi_yrs(:,:,end-recordYrs+1:end),3);
                Bw=nanmean(Bw_yrs(:,:,end-recordYrs+1:end,TempScenario),3);
                gainBw=nanmean(gainBw_yrs(:,:,end-recordYrs+1:end,TempScenario),3);
                Zw=nanmean(Zw_yrs(:,:,end-recordYrs+1:end,TempScenario),3);
                gainZw=nanmean(gainZw_yrs(:,:,end-recordYrs+1:end,TempScenario),3);
                TLallw=nanmean(TLallw_yrs(end-recordYrs+1:end,TempScenario));
                TLiw=nanmean(TLiw_yrs(:,:,end-recordYrs+1:end,TempScenario),3);
                
                %                 %Redo single-species LV hindcast fits by finding the best biomass
                %                 %and production quantiles to match model projections to
                %                 %forecast:
                %                 if sum(B(:))<eps
                %                     BLV1=nanmean(BLV4_yrs(:,:,end-recordYrs+1:end),3); %BLV, BLV1-6
                %                     gainBLV1=nanmean(gainBLV4_yrs(:,:,end-recordYrs+1:end),3);
                %                     BLV1w=nanmean(BLV4w_yrs(:,:,end-recordYrs+1:end),3);
                %                     gainBLV1w=nanmean(gainBLV4w_yrs(:,:,end-recordYrs+1:end),3); %gainBLV1w=nanmean(gainBLV4w_yrs(:,:,end-recordYrs+1:end,TempScenario),3);
                %                     BPquantiles=[BPquantiles;fitCode(4,:)];
                %                     BPpercDiff=[BPpercDiff;[NaN NaN]];
                %                 else
                %                     meanLV1B=mean(nansum(nanmean(BLV1_yrs(:,:,end-recordYrs+1:end),3),2));
                %                     meanLV2B=mean(nansum(nanmean(BLV2_yrs(:,:,end-recordYrs+1:end),3),2));
                %                     meanLV3B=mean(nansum(nanmean(BLV3_yrs(:,:,end-recordYrs+1:end),3),2));
                %                     meanLV4B=mean(nansum(nanmean(BLV4_yrs(:,:,end-recordYrs+1:end),3),2));
                %                     meanB=mean(nansum(B)); %mean observed biomass
                %                     meanLV1P=mean(nansum(nanmean(gainBLV1_yrs(:,:,end-recordYrs+1:end),3),2));
                %                     meanLV2P=mean(nansum(nanmean(gainBLV2_yrs(:,:,end-recordYrs+1:end),3),2));
                %                     meanLV3P=mean(nansum(nanmean(gainBLV3_yrs(:,:,end-recordYrs+1:end),3),2));
                %                     meanLV4P=mean(nansum(nanmean(gainBLV4_yrs(:,:,end-recordYrs+1:end),3),2));
                %                     meanP=mean(nansum(gainB)); %mean observed production
                %                     X=[ones(length(fitCode),1) fitCode fitCode(:,1).*fitCode(:,2)]; %set up regression predictor matrix
                %                     bB=regress([meanLV1B meanLV2B meanLV3B meanLV4B]',X); %regression model for B as function of B and P quantiles
                %                     bP=regress([meanLV1P meanLV2P meanLV3P meanLV4P]',X); %regression model for P as function of B and P quantiles
                %                     [fitCodeEsts,fval,exitflag,output] = fminsearchbnd(@(params) quantileFit(meanB,meanP,bB,bP,params),freeparstart,freeparmin,freeparmax); %fit quantiles to meanB and meanP
                %                     [percDiffB,percDiffP] = quantileOutputs(meanB,meanP,bB,bP,fitCodeEsts); %get percent differences between predicted and observed biomass and production
                %                     [r5,a5,z5,K5,flag5,raR25,r_T5,K_T5,K_T_ratio5,r_T_ratio5]=estSingleSpeciesModelmsy(Btrans,dBtrans,gainBtrans,P,fitCodeEsts); %fit growth model to all patches at once using estimated best quantiles
                %                     BPquantiles=[BPquantiles;fitCodeEsts];
                %                     BPpercDiff=[BPpercDiff;[percDiffB,percDiffP]];
                %                     disp(['case ' num2str(CaseNumber) ' iteration ' num2str(iteration) ': quantiles ' num2str(fitCodeEsts,2) ', %diff ' num2str([percDiffB,percDiffP],2)])
                %
                %                     %rerun single-species projection:
                %                     TimePts=[1:365:200*365+1]; %record every year
                %                     numPts=length(TimePts); %number of time points
                %                     gainBLV5=zeros(P.nx,P.n);
                %                     dBLV5=zeros(P.nx,P.n);
                %                     gainBLV5w=zeros(P.nx,P.n);
                %                     dBLV5w=zeros(P.nx,P.n);
                %                     BLV5w_yrs=zeros(P.nx,P.n,numPts);
                %                     gainBLV5w_yrs=zeros(P.nx,P.n,numPts);
                %                     BLV5=BLV1_yrs(:,:,1); %no warming case under estimated single-species dynamics (start at identical point as other simulations)
                %                     BLV5w=BLV5; %6C warming case under estimated single-species dynamics (start at identical point as other simulations)
                %                     for t = 1:TimePts(end) %run for 200 years (with daily time steps)
                %                         BLV5(BLV5<eps) = 0;
                %                         BLV5w(BLV5w<eps)= 0;
                %                         T1      = P.T;% + t.*P.dT; %<<< add this when time is right
                %                         T1w      = P.T + (t-1)*P.dT(3);
                %                         BLV5=sub_move(BLV5,P); %single species model no temp change, move
                %                         [gainBLV5 dBLV5] = sub_demogLVmsy(BLV5,T1,r5,a5,z5,P); % grow/die
                %                         BLV5w=sub_move(BLV5w,P); %single species model with temp change, move
                %                         [gainBLV5w dBLV5w] = sub_demogLVmsy(BLV5w,T1w,r5,a5,z5,P); % grow/die
                %                         tpos=find(t==TimePts);
                %                         if ~isempty(tpos)
                %                             BLV5_yrs(:,:,tpos)=BLV5;
                %                             gainBLV5_yrs(:,:,tpos)=gainBLV5;
                %                         end
                %                         BLV5 = BLV5 + dBLV5;
                %                         BLV5w = BLV5w + dBLV5w;
                %                     end
                
                %---------------Change assignments here for different single-species projections---------
                BLV1=nanmean(BLV4_yrs(:,:,end-recordYrs+1:end),3); %BLV, BLV1-6
                gainBLV1=nanmean(gainBLV4_yrs(:,:,end-recordYrs+1:end),3);
                BLV1w=nanmean(BLV4w_yrs(:,:,end-recordYrs+1:end),3);
                gainBLV1w=nanmean(gainBLV4w_yrs(:,:,end-recordYrs+1:end),3); %gainBLV1w=nanmean(gainBLV4w_yrs(:,:,end-recordYrs+1:end,TempScenario),3);
                %----------------------------------------------------------------------------------------
                %    end
                BtransEnd=nanmean(Btrans(:,:,end-recordYrs+1:end),3); %biomasses at the end of transcient period
                BtransEnd(BtransEnd<=eps)=nan;
                Bnan=B;
                Bnan(Bnan<=eps)=nan;
                BLV1nan=BLV1;
                BLV1nan(BLV1nan<=eps)=nan;
                Bwnan=Bw;
                Bwnan(Bwnan<=eps)=nan;
                BLV1wnan=BLV1w;
                BLV1wnan(BLV1wnan<=eps)=nan;
                
                numPatches=size(B,1);
                
                %variables for no warming/warming (col1, col2 unless otherwise commented) cases:
                %global (for heterotrophs)
                %totBiomass=[totBiomass; [nansum(B(:)) nansum(Bw(:)) nansum(BLV1(:)) nansum(BLV1w(:))]]; %total biomass
                totBiomass=[totBiomass; [mean(nansum(B,2)) mean(nansum(Bw,2)) mean(nansum(BLV1,2)) mean(nansum(BLV1w,2))]]; %total biomass (as mean biomass/m^3)
                allBody_Biomass=cat(3,allBody_Biomass,[P.S(2:end);nansum(B);nansum(Bw);nansum(BLV1);nansum(BLV1w)]);
                totProd=[totProd; [mean(nansum(gainB,2)) mean(nansum(gainBw,2)) mean(nansum(gainBLV1,2)) mean(nansum(gainBLV1w,2))]]; %total production (as mean production/(day m^3))
                %totRich=[totRich; [nansum(nansum(B)>eps*numPatches) nansum(nansum(Bw)>eps*numPatches) nansum(nansum(BLV1)>eps*numPatches) nansum(nansum(BLV1w)>eps*numPatches)]]; %species richness
                %totBeta=[totBeta; [nansum(nansum(B)>eps*numPatches)/nanmean(nansum(B>eps*numPatches,2)) nansum(nansum(Bw)>eps*numPatches)/nanmean(nansum(Bw>eps*numPatches,2)) nansum(nansum(BLV1)>eps*numPatches)/nanmean(nansum(BLV1>eps*numPatches,2)) nansum(nansum(BLV1w)>eps*numPatches)/nanmean(nansum(BLV1w>eps*numPatches,2))]]; %beta diversity
                %for LV single-species model projection, count species as
                %extinct if biomass<eps*numPatches*ExtinctMag
                ExtinctMag=1;
                totRich=[totRich; [nansum(nansum(B)>eps*numPatches) nansum(nansum(Bw)>eps*numPatches) nansum(nansum(BLV1)>eps*numPatches*ExtinctMag) nansum(nansum(BLV1w)>eps*numPatches*ExtinctMag)]]; %species richness
                totBeta=[totBeta; [nansum(nansum(B)>eps*numPatches)/nanmean(nansum(B>eps*numPatches,2)) nansum(nansum(Bw)>eps*numPatches)/nanmean(nansum(Bw>eps*numPatches,2)) nansum(nansum(BLV1)>eps*numPatches*ExtinctMag)/nanmean(nansum(BLV1>eps*numPatches*ExtinctMag,2)) nansum(nansum(BLV1w)>eps*numPatches*ExtinctMag)/nanmean(nansum(BLV1w>eps*numPatches*ExtinctMag,2))]]; %beta diversity
                totAlpha=[totAlpha; [nanmean(nansum(B>eps*numPatches,2)) nanmean(nansum(Bw>eps*numPatches,2)) nanmean(nansum(BLV1>eps*numPatches*ExtinctMag,2)) nanmean(nansum(BLV1w>eps*numPatches*ExtinctMag,2))]]; %alpha diversity
                totTrophicLevel=[totTrophicLevel; [TLall TLallw]]; %mean global trophic level
                maxTrophicLevel=[maxTrophicLevel; [max(TLi) max(TLiw)]]; %maximum non-extinct global trophic species
                totBodyMass=[totBodyMass; [nansum(nansum(B).*P.S(2:end))/nansum(B(:)) nansum(nansum(Bw).*P.S(2:end))/nansum(Bw(:)) nansum(nansum(BLV1).*P.S(2:end))/nansum(BLV1(:)) nansum(nansum(BLV1w).*P.S(2:end))/nansum(BLV1w(:))]]; %mean body mass (by biomass)
                tempMaxBodyMass=P.S(max(find(nansum(B)>eps*numPatches))+1);
                if isempty(tempMaxBodyMass)
                    tempMaxBodyMass=NaN;
                end
                tempMaxBodyMassW=P.S(max(find(nansum(Bw)>eps*numPatches))+1);
                if isempty(tempMaxBodyMassW)
                    tempMaxBodyMassW=NaN;
                end
                tempMaxBodyMassLV1=P.S(max(find(nansum(BLV1)>0))+1);
                if isempty(tempMaxBodyMassLV1)
                    tempMaxBodyMassLV1=NaN;
                end
                tempMaxBodyMassLV1w=P.S(max(find(nansum(BLV1w)>0))+1);
                if isempty(tempMaxBodyMassLV1w)
                    tempMaxBodyMassLV1w=NaN;
                end
                maxBodyMass=[maxBodyMass; [tempMaxBodyMass tempMaxBodyMassW tempMaxBodyMassLV1 tempMaxBodyMassLV1w]]; %maximum body mass
                
                [maxB findMaxB]=max(nansum(B));
                [maxBLV1 findMaxBLV1]=max(nansum(BLV1));
                [maxP findMaxP]=max(nansum(gainB));
                [maxPLV1 findMaxPLV1]=max(nansum(gainBLV1));
                [maxBw findMaxBw]=max(nansum(Bw));
                [maxBLV1w findMaxBLV1w]=max(nansum(BLV1w));
                [maxPw findMaxPw]=max(nansum(gainBw));
                [maxPLV1w findMaxPLV1w]=max(nansum(gainBLV1w));
                maxBiomass=[maxBiomass; [maxB./nansum(B(:)) maxBw./nansum(Bw(:)) maxBLV1./nansum(BLV1(:)) maxBLV1w./nansum(BLV1w(:))]]; %fraction
                maxProd=[maxProd; [maxP./nansum(gainB(:)) maxPw./nansum(gainBw(:)) maxPLV1./nansum(gainBLV1(:)) maxPLV1w./nansum(gainBLV1w(:))]]; %fraction (100% if total community production<0)
                maxProd(maxProd<0)=1;
                %maxProd=[maxProd; [maxP maxPw maxPLV1 maxPLV1w]];
                maxBiomassBS=[maxBiomassBS; [P.S(findMaxB+1) P.S(findMaxBw+1) P.S(findMaxBLV1+1) P.S(findMaxBLV1w+1)]];
                maxProdBS=[maxProdBS; [P.S(findMaxP+1) P.S(findMaxPw+1) P.S(findMaxPLV1+1) P.S(findMaxPLV1w+1)]];
                consResRatio=[consResRatio; nansum(B(:))/nansum(Z(:)) nansum(Bw(:))/nansum(Zw(:))];
                fractSpeciesProd=[fractSpeciesProd; [nansum(nansum(gainB)>0)/nansum(nansum(B)>eps*numPatches) nansum(nansum(gainBw)>0)/nansum(nansum(Bw)>eps*numPatches) nansum(nansum(gainBLV1)>0)/nansum(nansum(BLV1)>eps*numPatches) nansum(nansum(gainBLV1w)>0)/nansum(nansum(BLV1w)>eps*numPatches)]];
                %range shift (median)
                [dummy BtransEndLoc]=max(BtransEnd);
                BtransEndLoc(isnan(dummy))=nan;
                originLoc=[originLoc; BtransEndLoc]; %location with highest biomass of each species at the end of transcient period
                [dummy BLoc]=max(Bnan);
                BLoc(isnan(dummy))=nan;
                finalLoc=[finalLoc; BLoc]; %location with highest biomass of each species at the end of warming period
                [dummy BLV1Loc]=max(BLV1nan);
                BLV1Loc(isnan(dummy))=nan;
                finalLocLV1=[finalLocLV1; BLV1Loc]; %location with highest biomass of each species at the end of warming period
                rangeShift0=[rangeShift0; nanmean(BLoc-BtransEndLoc) nanmean(BLV1Loc-BtransEndLoc)]; %average location shift
                [dummy originBiomLoc]=max(nansum(BtransEnd,2));
                [dummy finalBiomLoc]=max(nansum(Bnan,2));
                [dummy finalBiomLV1Loc]=max(nansum(BLV1nan,2));
                BiomShift0=[BiomShift0; finalBiomLoc-originBiomLoc finalBiomLV1Loc-originBiomLoc];
                
                [dummy BwLoc]=max(Bwnan);
                BwLoc(isnan(dummy))=nan;
                finalwLoc=[finalwLoc; BwLoc]; %location with highest biomass of each species at the end of warming period
                [dummy BLV1wLoc]=max(BLV1wnan);
                BLV1wLoc(isnan(dummy))=nan;
                finalwLocLV1=[finalwLocLV1; BLV1wLoc]; %location with highest biomass of each species at the end of warming period
                rangeShift=[rangeShift; nanmean(BwLoc-BtransEndLoc) nanmean(BLV1wLoc-BtransEndLoc)]; %average location shift
                [dummy finalwBiomLoc]=max(nansum(Bwnan,2));
                [dummy finalwBiomLV1Loc]=max(nansum(BLV1wnan,2));
                BiomShift=[BiomShift; finalwBiomLoc-originBiomLoc finalwBiomLV1Loc-originBiomLoc];
                All_Shifts=cat(3,All_Shifts,[BLoc-BtransEndLoc;BwLoc-BtransEndLoc;BLV1Loc-BtransEndLoc;BLV1wLoc-BtransEndLoc]);
                
                %range shift (centroid or mean location)
                meanInitCentroid=nansum(BtransEnd.*[1:11]')./nansum(BtransEnd);
                meanInitCentroid(nansum(BtransEnd)==0)=NaN;
                InitCentroid=[InitCentroid; meanInitCentroid]; %mean location of each species at the end of transcient period
                meanFinalCentroid=nansum(Bnan.*[1:11]')./nansum(Bnan);
                meanFinalCentroid(nansum(Bnan)==0)=NaN;
                meanFinalLVCentroid=nansum(BLV1nan.*[1:11]')./nansum(BLV1nan);
                meanFinalLVCentroid(nansum(BLV1nan)==0)=NaN;
                FinalCentroid=[FinalCentroid; meanFinalCentroid]; %mean location of each species at the end of warming period
                FinalLVCentroid=[FinalLVCentroid; meanFinalLVCentroid]; %projected mean location of each species at the end of warming period
                meanFinalwCentroid=nansum(Bwnan.*[1:11]')./nansum(Bwnan);
                meanFinalwCentroid(nansum(Bwnan)==0)=NaN;
                meanFinalwLVCentroid=nansum(BLV1wnan.*[1:11]')./nansum(BLV1wnan);
                meanFinalwLVCentroid(nansum(BLV1wnan)==0)=NaN;
                FinalwCentroid=[FinalwCentroid; meanFinalwCentroid]; %mean location of each species at the end of warming period
                FinalwLVCentroid=[FinalwLVCentroid; meanFinalwLVCentroid]; %projected mean location of each species at the end of warming period
                CentroidShift0=[CentroidShift0; nanmean(meanFinalCentroid-meanInitCentroid) nanmean(meanFinalLVCentroid-meanInitCentroid)]; %average centroid shift
                CentroidShift=[CentroidShift; nanmean(meanFinalwCentroid-meanInitCentroid) nanmean(meanFinalwLVCentroid-meanInitCentroid)]; %average centroid shift
                All_CentroidShifts=cat(3,All_CentroidShifts,[meanFinalCentroid-meanInitCentroid;meanFinalwCentroid-meanInitCentroid;meanFinalLVCentroid-meanInitCentroid;meanFinalwLVCentroid-meanInitCentroid]);
                
                %trailing edge shift (bootstrap from sample distribution ~sample size=1000*densities)
                trailingQ=0.025; %lower quantile
                leadingQ=0.975; %upper quantile
                minSpeciesBiomass=trailingQ*nansum(BtransEnd); %determine trailing biomass density
                maxSpeciesBiomass=leadingQ*nansum(BtransEnd); %determine leading biomass density
                tempInitTrailing=sum(cumsum(BtransEnd,'omitnan')<minSpeciesBiomass)+1;
                tempInitTrailing(nansum(BtransEnd)==0)=NaN;
                tempInitLeading=sum(cumsum(BtransEnd,'omitnan')<maxSpeciesBiomass)+1;
                tempInitLeading(nansum(BtransEnd)==0)=NaN;
                InitTrailing=[InitTrailing; tempInitTrailing]; %mean location of each species at the end of transcient period
                InitLeading=[InitLeading; tempInitLeading]; %mean location of each species at the end of transcient period
                
                minSpeciesBiomass=trailingQ*nansum(Bnan); %determine trailing biomass density
                maxSpeciesBiomass=leadingQ*nansum(Bnan); %determine leading biomass density
                tempFinalTrailing=sum(cumsum(Bnan,'omitnan')<minSpeciesBiomass)+1;
                tempFinalTrailing(nansum(Bnan)==0)=NaN;
                tempFinalLeading=sum(cumsum(Bnan,'omitnan')<maxSpeciesBiomass)+1;
                tempFinalLeading(nansum(Bnan)==0)=NaN;
                FinalTrailing=[FinalTrailing; tempFinalTrailing]; %mean location of each species at the end of transcient period
                FinalLeading=[FinalLeading; tempFinalLeading]; %mean location of each species at the end of transcient period
                
                minSpeciesBiomass=trailingQ*nansum(Bwnan); %determine trailing biomass density
                maxSpeciesBiomass=leadingQ*nansum(Bwnan); %determine leading biomass density
                tempFinalwTrailing=sum(cumsum(Bwnan,'omitnan')<minSpeciesBiomass)+1;
                tempFinalwTrailing(nansum(Bwnan)==0)=NaN;
                tempFinalwLeading=sum(cumsum(Bwnan,'omitnan')<maxSpeciesBiomass)+1;
                tempFinalwLeading(nansum(Bwnan)==0)=NaN;
                FinalwTrailing=[FinalwTrailing; tempFinalwTrailing]; %mean location of each species at the end of transcient period
                FinalwLeading=[FinalwLeading; tempFinalwLeading]; %mean location of each species at the end of transcient period
                
                minSpeciesBiomass=trailingQ*nansum(BLV1nan); %determine trailing biomass density
                maxSpeciesBiomass=leadingQ*nansum(BLV1nan); %determine leading biomass density
                tempFinalLVTrailing=sum(cumsum(BLV1nan,'omitnan')<minSpeciesBiomass)+1;
                tempFinalLVTrailing(nansum(BLV1nan)==0)=NaN;
                tempFinalLVLeading=sum(cumsum(BLV1nan,'omitnan')<maxSpeciesBiomass)+1;
                tempFinalLVLeading(nansum(BLV1nan)==0)=NaN;
                FinalLVTrailing=[FinalLVTrailing; tempFinalLVTrailing]; %mean location of each species at the end of transcient period
                FinalLVLeading=[FinalLVLeading; tempFinalLVLeading]; %mean location of each species at the end of transcient period
                
                minSpeciesBiomass=trailingQ*nansum(BLV1wnan); %determine trailing biomass density
                maxSpeciesBiomass=leadingQ*nansum(BLV1wnan); %determine leading biomass density
                tempFinalwLVTrailing=sum(cumsum(BLV1wnan,'omitnan')<minSpeciesBiomass)+1;
                tempFinalwLVTrailing(nansum(BLV1wnan)==0)=NaN;
                tempFinalwLVLeading=sum(cumsum(BLV1wnan,'omitnan')<maxSpeciesBiomass)+1;
                tempFinalwLVLeading(nansum(BLV1wnan)==0)=NaN;
                FinalwLVTrailing=[FinalwLVTrailing; tempFinalwLVTrailing]; %mean location of each species at the end of transcient period
                FinalwLVLeading=[FinalwLVLeading; tempFinalwLVLeading]; %mean location of each species at the end of transcient period
                
                TrailingShift0=[TrailingShift0; nanmean(tempFinalTrailing-tempInitTrailing) nanmean(tempFinalLVTrailing-tempInitTrailing)]; %average trailing edge shift under no warming
                TrailingShift=[TrailingShift; nanmean(tempFinalwTrailing-tempInitTrailing) nanmean(tempFinalwLVTrailing-tempInitTrailing)]; %average trailing edge shift under warming
                LeadingShift0=[LeadingShift0; nanmean(tempFinalLeading-tempInitLeading) nanmean(tempFinalLVLeading-tempInitLeading)]; %average trailing edge shift under no warming
                LeadingShift=[LeadingShift; nanmean(tempFinalwLeading-tempInitLeading) nanmean(tempFinalwLVLeading-tempInitLeading)]; %average trailing edge shift under warming
                All_TrailingShifts=cat(3,All_TrailingShifts,[tempFinalTrailing-tempInitTrailing;tempFinalwTrailing-tempInitTrailing;tempFinalLVTrailing-tempInitTrailing;tempFinalwLVTrailing-tempInitTrailing]);
                All_LeadingShifts=cat(3,All_LeadingShifts,[tempFinalLeading-tempInitLeading;tempFinalwLeading-tempInitLeading;tempFinalLVLeading-tempInitLeading;tempFinalwLVLeading-tempInitLeading]);
                
            end
        end
    end
end



%numIt=length(paramIndices)/numCases; %number of iterations per case
%(assuming equal numbers for each case)
tempChanges=P.dT*73000; %temperature changes corresponding to temperature scenario
dispersalLabels={'0' '1' '10^3' '10^6' '10^9' '10^{12}'};
scrsz = get(0,'ScreenSize');

FileName=sprintf('WarmingMovementStats_%s.mat', TimeData);
save(FileName);

%figs(1)=figure ('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/2 scrsz(4)]);
%subplot(4,2,1)
figs(1)=figure ('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/2 scrsz(4)/4]);
subplot(1,2,1)
hold on
bl1=boundedline([1:numCases], nanmean(reshape((totBiomass(:,1)),numIt,[])),[nanstd(reshape((totBiomass(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((totBiomass(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2);
%bl1=boundedline([1:numCases], nanmean(reshape((totBiomass(:,2)),numIt,[])),[nanstd(reshape((totBiomass(:,2)),numIt,[]))*1.96./((sum(~isnan(reshape((totBiomass(:,2)),numIt,[])))).^0.5)]','r','alpha'); drawnow; set(bl1,'linewidth',2);
bl1=boundedline([1:numCases], nanmean(reshape((totBiomass(:,3)),numIt,[])),[nanstd(reshape((totBiomass(:,3)),numIt,[]))*1.96./((sum(~isnan(reshape((totBiomass(:,3)),numIt,[])))).^0.5)]','--k','alpha','transparency', 0.1);
%bl1=boundedline([1:numCases], nanmean(reshape((totBiomass(:,4)),numIt,[])),[nanstd(reshape((totBiomass(:,4)),numIt,[]))*1.96./((sum(~isnan(reshape((totBiomass(:,4)),numIt,[])))).^0.5)]','--r','alpha','transparency', 0.1);
%bl1=boundedline([1:numCases], nanmean(reshape(log10(totBiomass(:,1)),numIt,[])),[nanstd(reshape(log10(totBiomass(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape(log10(totBiomass(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2);
%bl1=boundedline([1:numCases], nanmean(reshape(log10(totBiomass(:,2)),numIt,[])),[nanstd(reshape(log10(totBiomass(:,2)),numIt,[]))*1.96./((sum(~isnan(reshape(log10(totBiomass(:,2)),numIt,[])))).^0.5)]','r','alpha'); drawnow; set(bl1,'linewidth',2);
%bl1=boundedline([1:numCases], nanmean(reshape(log10(totBiomass(:,3)),numIt,[])),[nanstd(reshape(log10(totBiomass(:,3)),numIt,[]))*1.96./((sum(~isnan(reshape(log10(totBiomass(:,3)),numIt,[])))).^0.5)]','--b','alpha','transparency', 0.1);
%bl1=boundedline([1:numCases], nanmean(reshape(log10(totBiomass(:,4)),numIt,[])),[nanstd(reshape(log10(totBiomass(:,4)),numIt,[]))*1.96./((sum(~isnan(reshape(log10(totBiomass(:,4)),numIt,[])))).^0.5)]','--r','alpha','transparency', 0.1);
xlim([1 numCases]); xticks([1:numCases]); xticklabels(dispersalLabels);
%title 'outcomes'
xlabel 'dispersal rate'
ylabel 'biomass [gm^{-3}]'
legend off
%subplot(4,2,2)
subplot(1,2,2)
hold on
%gscatter([paramIndices;paramIndices],totProd(:),warmingIndices(:),'br','o')
hold on
bl1=boundedline([1:numCases], nanmean(reshape((totProd(:,1)),numIt,[])),[nanstd(reshape((totProd(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((totProd(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2);
%bl1=boundedline([1:numCases], nanmean(reshape((totProd(:,2)),numIt,[])),[nanstd(reshape((totProd(:,2)),numIt,[]))*1.96./((sum(~isnan(reshape((totProd(:,2)),numIt,[])))).^0.5)]','r','alpha'); drawnow; set(bl1,'linewidth',2);
bl1=boundedline([1:numCases], nanmean(reshape((totProd(:,3)),numIt,[])),[nanstd(reshape((totProd(:,3)),numIt,[]))*1.96./((sum(~isnan(reshape((totProd(:,3)),numIt,[])))).^0.5)]','--k','alpha','transparency', 0.1);
%bl1=boundedline([1:numCases], nanmean(reshape((totProd(:,4)),numIt,[])),[nanstd(reshape((totProd(:,4)),numIt,[]))*1.96./((sum(~isnan(reshape((totProd(:,4)),numIt,[])))).^0.5)]','--r','alpha','transparency', 0.1);
xlim([1 numCases]); xticks([1:numCases]); xticklabels(dispersalLabels);
%title (['warming ' num2str(tempChanges(TempScenario)) '\circC'])
xlabel 'dispersal rate'
ylabel 'production [gm^-{3}/day]'
legend off
% subplot(4,2,3)
% %gscatter([paramIndices;paramIndices],totRich(:),warmingIndices(:),'br','o')
% hold on
% bl1=boundedline([1:numCases], nanmean(reshape((totRich(:,1)),numIt,[])),[nanstd(reshape((totRich(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((totRich(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2);
% bl1=boundedline([1:numCases], nanmean(reshape((totRich(:,2)),numIt,[])),[nanstd(reshape((totRich(:,2)),numIt,[]))*1.96./((sum(~isnan(reshape((totRich(:,2)),numIt,[])))).^0.5)]','r','alpha'); drawnow; set(bl1,'linewidth',2);
% %bl1=boundedline([1:numCases], nanmean(reshape((totRich(:,3)),numIt,[])),[nanstd(reshape((totRich(:,3)),numIt,[]))*1.96./((sum(~isnan(reshape((totRich(:,3)),numIt,[])))).^0.5)]','--b','alpha','transparency', 0.1);
% %bl1=boundedline([1:numCases], nanmean(reshape((totRich(:,4)),numIt,[])),[nanstd(reshape((totRich(:,4)),numIt,[]))*1.96./((sum(~isnan(reshape((totRich(:,4)),numIt,[])))).^0.5)]','--r','alpha','transparency', 0.1);
% xlim([1 numCases]); xticks([1:numCases]); xticklabels(dispersalLabels);
% xlabel 'dispersal rate'
% ylabel 'richness'
% title ''
% legend off
% subplot(4,2,4)
% %gscatter([paramIndices;paramIndices],totBeta(:),warmingIndices(:),'br','o')
% hold on
% %bl1=boundedline([1:numCases], nanmean(reshape((totAlpha(:,1)),numIt,[])),[nanstd(reshape((totAlpha(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((totAlpha(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2);
% %bl1=boundedline([1:numCases], nanmean(reshape((totAlpha(:,2)),numIt,[])),[nanstd(reshape((totAlpha(:,2)),numIt,[]))*1.96./((sum(~isnan(reshape((totAlpha(:,2)),numIt,[])))).^0.5)]','r','alpha'); drawnow; set(bl1,'linewidth',2);
% bl1=boundedline([1:numCases], nanmean(reshape((totBeta(:,1)),numIt,[])),[nanstd(reshape((totBeta(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((totBeta(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2);
% bl1=boundedline([1:numCases], nanmean(reshape((totBeta(:,2)),numIt,[])),[nanstd(reshape((totBeta(:,2)),numIt,[]))*1.96./((sum(~isnan(reshape((totBeta(:,2)),numIt,[])))).^0.5)]','r','alpha'); drawnow; set(bl1,'linewidth',2);
% %bl1=boundedline([1:numCases], nanmean(reshape((totBeta(:,3)),numIt,[])),[nanstd(reshape((totBeta(:,3)),numIt,[]))*1.96./((sum(~isnan(reshape((totBeta(:,3)),numIt,[])))).^0.5)]','--b','alpha','transparency', 0.1);
% %bl1=boundedline([1:numCases], nanmean(reshape((totBeta(:,4)),numIt,[])),[nanstd(reshape((totBeta(:,4)),numIt,[]))*1.96./((sum(~isnan(reshape((totBeta(:,4)),numIt,[])))).^0.5)]','--r','alpha','transparency', 0.1);
% xlim([1 numCases]); xticks([1:numCases]); xticklabels(dispersalLabels);
% xlabel 'dispersal rate'
% ylabel 'beta diversity'
% %ylabel 'mean alpha diversity'
% title ''
% legend off
% subplot(4,2,7)
% %gscatter([paramIndices;paramIndices],totTrophicLevel(:),warmingIndices(:),'br','o')
% hold on
% bl1=boundedline([1:numCases], nanmean(reshape((totTrophicLevel(:,1)),numIt,[])),[nanstd(reshape((totTrophicLevel(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((totTrophicLevel(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2);
% bl1=boundedline([1:numCases], nanmean(reshape((totTrophicLevel(:,2)),numIt,[])),[nanstd(reshape((totTrophicLevel(:,2)),numIt,[]))*1.96./((sum(~isnan(reshape((totTrophicLevel(:,2)),numIt,[])))).^0.5)]','r','alpha'); drawnow; set(bl1,'linewidth',2);
% xlim([1 numCases]); xticks([1:numCases]); xticklabels(dispersalLabels);
% xlabel 'dispersal rate'
% ylabel 'mean trophic level'
% title ''
% legend off
% subplot(4,2,8)
% %gscatter([paramIndices;paramIndices],totTrophicLevel(:),warmingIndices(:),'br','o')
% hold on
% bl1=boundedline([1:numCases], nanmean(reshape((maxTrophicLevel(:,1)),numIt,[])),[nanstd(reshape((maxTrophicLevel(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((totTrophicLevel(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2);
% bl1=boundedline([1:numCases], nanmean(reshape((maxTrophicLevel(:,2)),numIt,[])),[nanstd(reshape((maxTrophicLevel(:,2)),numIt,[]))*1.96./((sum(~isnan(reshape((totTrophicLevel(:,2)),numIt,[])))).^0.5)]','r','alpha'); drawnow; set(bl1,'linewidth',2);
% xlim([1 numCases]); xticks([1:numCases]); xticklabels(dispersalLabels);
% xlabel 'dispersal rate'
% ylabel 'max trophic level'
% title ''
% legend off
% subplot(4,2,5)
% %gscatter([paramIndices;paramIndices],totTrophicLevel(:),warmingIndices(:),'br','o')
% hold on
% bl1=boundedline([1:numCases], nanmean(reshape((totBodyMass(:,1)),numIt,[])),[nanstd(reshape((totBodyMass(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((totBodyMass(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2);
% bl1=boundedline([1:numCases], nanmean(reshape((totBodyMass(:,2)),numIt,[])),[nanstd(reshape((totBodyMass(:,2)),numIt,[]))*1.96./((sum(~isnan(reshape((totBodyMass(:,2)),numIt,[])))).^0.5)]','r','alpha'); drawnow; set(bl1,'linewidth',2);
% %bl1=boundedline([1:numCases], nanmean(reshape((totBodyMass(:,3)),numIt,[])),[nanstd(reshape((totBodyMass(:,3)),numIt,[]))*1.96./((sum(~isnan(reshape((totBodyMass(:,3)),numIt,[])))).^0.5)]','--b','alpha','transparency', 0.1);
% %bl1=boundedline([1:numCases], nanmean(reshape((totBodyMass(:,4)),numIt,[])),[nanstd(reshape((totBodyMass(:,4)),numIt,[]))*1.96./((sum(~isnan(reshape((totBodyMass(:,4)),numIt,[])))).^0.5)]','--r','alpha','transparency', 0.1);
% xlim([1 numCases]); xticks([1:numCases]); xticklabels(dispersalLabels);
% xlabel 'dispersal rate'
% ylabel 'mean log_{10}(body size [g])'
% title ''
% legend off
% subplot(4,2,6)
% %gscatter([paramIndices;paramIndices],totTrophicLevel(:),warmingIndices(:),'br','o')
% hold on
% % bl1=boundedline([1:numCases], nanmean(reshape((maxBodyMass(:,1)),numIt,[])),[nanstd(reshape((maxBodyMass(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((maxBodyMass(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2);
% % bl1=boundedline([1:numCases], nanmean(reshape((maxBodyMass(:,2)),numIt,[])),[nanstd(reshape((maxBodyMass(:,2)),numIt,[]))*1.96./((sum(~isnan(reshape((maxBodyMass(:,2)),numIt,[])))).^0.5)]','r','alpha'); drawnow; set(bl1,'linewidth',2);
% % bl1=boundedline([1:numCases], nanmean(reshape((maxBodyMass(:,3)),numIt,[])),[nanstd(reshape((maxBodyMass(:,3)),numIt,[]))*1.96./((sum(~isnan(reshape((maxBodyMass(:,3)),numIt,[])))).^0.5)]','--b','alpha','transparency', 0.1);
% % bl1=boundedline([1:numCases], nanmean(reshape((maxBodyMass(:,4)),numIt,[])),[nanstd(reshape((maxBodyMass(:,4)),numIt,[]))*1.96./((sum(~isnan(reshape((maxBodyMass(:,4)),numIt,[])))).^0.5)]','--r','alpha','transparency', 0.1);
% bl1=boundedline([1:numCases], nanmean(reshape((consResRatio(:,1)),numIt,[])),[nanstd(reshape((consResRatio(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((consResRatio(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2);
% bl1=boundedline([1:numCases], nanmean(reshape((consResRatio(:,2)),numIt,[])),[nanstd(reshape((consResRatio(:,2)),numIt,[]))*1.96./((sum(~isnan(reshape((consResRatio(:,2)),numIt,[])))).^0.5)]','r','alpha'); drawnow; set(bl1,'linewidth',2);
% xlim([1 numCases]); xticks([1:numCases]); xticklabels(dispersalLabels);
% xlabel 'dispersal rate'
% %ylabel 'max body size'
% ylabel 'consumer:resource'
% title ''
% legend off
%
% figs(2)=figure ('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/2 scrsz(4)]);
% subplot(4,2,1)
% hold on
% %violin(reshape(log2(totBiomass(:,2)./totBiomass(:,1)),numIt,[]),'facecolor','w','facealpha',1,'mc',[],'medc',[]);
% bl1=boundedline([1:numCases], nanmean(reshape((totBiomass(:,2)-totBiomass(:,1)),numIt,[])),[nanstd(reshape((totBiomass(:,2)-totBiomass(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((totBiomass(:,2)-totBiomass(:,1)),numIt,[])))).^0.5)]','r','alpha'); drawnow; set(bl1,'linewidth',2);
% bl1=boundedline([1:numCases], nanmean(reshape((totBiomass(:,4)-totBiomass(:,3)),numIt,[])),[nanstd(reshape((totBiomass(:,4)-totBiomass(:,3)),numIt,[]))*1.96./((sum(~isnan(reshape((totBiomass(:,4)-totBiomass(:,3)),numIt,[])))).^0.5)]','--k','alpha');
% refl=refline(0,0);
% set(refl,'color','k')
% xlim([1 numCases]); xticks([1:numCases]); xticklabels(dispersalLabels);
% title 'warming effects'
% xlabel 'dispersal rate'
% ylabel '\Delta biomass [gm^{-3}]'
% subplot(4,2,2)
% %violin(reshape((totProd(:,2)-totProd(:,1)),numIt,[]),'facecolor','w','facealpha',1,'mc',[],'medc',[]);
% hold on
% bl1=boundedline([1:numCases], nanmean(reshape((totProd(:,2)-totProd(:,1)),numIt,[])),[nanstd(reshape((totProd(:,2)-totProd(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((totProd(:,2)-totProd(:,1)),numIt,[])))).^0.5)]','r','alpha'); drawnow; set(bl1,'linewidth',2);
% bl1=boundedline([1:numCases], nanmean(reshape((totProd(:,4)-totProd(:,3)),numIt,[])),[nanstd(reshape((totProd(:,4)-totProd(:,3)),numIt,[]))*1.96./((sum(~isnan(reshape((totProd(:,4)-totProd(:,3)),numIt,[])))).^0.5)]','--k','alpha');
% refl=refline(0,0);
% set(refl,'color','k')
% xlim([1 numCases]); xticks([1:numCases]); xticklabels(dispersalLabels);
% title (['warming ' num2str(tempChanges(TempScenario)) '\circC'])
% xlabel 'dispersal rate'
% ylabel '\Delta production [gm^{-3}/day]'
%
% subplot(4,2,3)
% %violin(reshape(log2(totRich(:,2)./totRich(:,1)),numIt,[]),'facecolor','w','facealpha',1,'mc',[],'medc',[]);
% hold on
% bl1=boundedline([1:numCases], nanmean(reshape((totRich(:,2)-totRich(:,1)),numIt,[])),[nanstd(reshape((totRich(:,2)-totRich(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((totRich(:,2)-totRich(:,1)),numIt,[])))).^0.5)]','r','alpha'); drawnow; set(bl1,'linewidth',2);
% bl1=boundedline([1:numCases], nanmean(reshape((totRich(:,4)-totRich(:,3)),numIt,[])),[nanstd(reshape((totRich(:,4)-totRich(:,3)),numIt,[]))*1.96./((sum(~isnan(reshape((totRich(:,4)-totRich(:,3)),numIt,[])))).^0.5)]','--k','alpha');
% refl=refline(0,0);
% set(refl,'color','k')
% xlim([1 numCases]); xticks([1:numCases]); xticklabels(dispersalLabels);
% xlabel 'dispersal rate'
% ylabel '\Delta richness'
% subplot(4,2,4)
% %violin(reshape(log2(totBeta(:,2)./totBeta(:,1)),numIt,[]),'facecolor','w','facealpha',1,'mc',[],'medc',[]);
% hold on
% %bl1=boundedline([1:numCases], nanmean(reshape((totAlpha(:,2)-totAlpha(:,1)),numIt,[])),[nanstd(reshape((totAlpha(:,2)-totAlpha(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((totAlpha(:,2)-totAlpha(:,1)),numIt,[])))).^0.5)]','r','alpha'); drawnow; set(bl1,'linewidth',2);
% %bl1=boundedline([1:numCases], nanmean(reshape((totAlpha(:,4)-totAlpha(:,3)),numIt,[])),[nanstd(reshape((totAlpha(:,4)-totAlpha(:,3)),numIt,[]))*1.96./((sum(~isnan(reshape((totAlpha(:,4)-totAlpha(:,3)),numIt,[])))).^0.5)]','--k','alpha');
% bl1=boundedline([1:numCases], nanmean(reshape((totBeta(:,2)-totBeta(:,1)),numIt,[])),[nanstd(reshape((totBeta(:,2)-totBeta(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((totBeta(:,2)-totBeta(:,1)),numIt,[])))).^0.5)]','r','alpha'); drawnow; set(bl1,'linewidth',2);
% bl1=boundedline([1:numCases], nanmean(reshape((totBeta(:,4)-totBeta(:,3)),numIt,[])),[nanstd(reshape((totBeta(:,4)-totBeta(:,3)),numIt,[]))*1.96./((sum(~isnan(reshape((totBeta(:,4)-totBeta(:,3)),numIt,[])))).^0.5)]','--k','alpha');
% refl=refline(0,0);
% set(refl,'color','k')
% xlim([1 numCases]); xticks([1:numCases]); xticklabels(dispersalLabels);
% xlabel 'dispersal rate'
% %ylabel '\Delta mean alpha diversity'
% ylabel '\Delta beta diversity'
% subplot(4,2,7)
% %violin(reshape(log2(totTrophicLevel(:,2)./totTrophicLevel(:,1)),numIt,[]),'facecolor','w','facealpha',1,'mc',[],'medc',[]);
% hold on
% bl1=boundedline([1:numCases], nanmean(reshape((totTrophicLevel(:,2)-totTrophicLevel(:,1)),numIt,[])),[nanstd(reshape((totTrophicLevel(:,2)-totTrophicLevel(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((totTrophicLevel(:,2)-totTrophicLevel(:,1)),numIt,[])))).^0.5)]','r','alpha'); drawnow; set(bl1,'linewidth',2);
% refl=refline(0,0);
% set(refl,'color','k')
% xlim([1 numCases]); xticks([1:numCases]); xticklabels(dispersalLabels);
% xlabel 'dispersal rate'
% ylabel '\Delta mean trophic level'
% subplot(4,2,8)
% %violin(reshape(log2(maxTrophicLevel(:,2)./maxTrophicLevel(:,1)),numIt,[]),'facecolor','w','facealpha',1,'mc',[],'medc',[]);
% hold on
% bl1=boundedline([1:numCases], nanmean(reshape((maxTrophicLevel(:,2)-maxTrophicLevel(:,1)),numIt,[])),[nanstd(reshape((maxTrophicLevel(:,2)-maxTrophicLevel(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((maxTrophicLevel(:,2)-maxTrophicLevel(:,1)),numIt,[])))).^0.5)]','r','alpha'); drawnow; set(bl1,'linewidth',2);
% refl=refline(0,0);
% set(refl,'color','k')
% xlim([1 numCases]); xticks([1:numCases]); xticklabels(dispersalLabels);
% xlabel 'dispersal rate'
% ylabel '\Delta max trophic level'
% subplot(4,2,5)
% %violin(reshape(log2(totBodyMass(:,2)./totBodyMass(:,1)),numIt,[]),'facecolor','w','facealpha',1,'mc',[],'medc',[]);
% hold on
% bl1=boundedline([1:numCases], nanmean(reshape((totBodyMass(:,2)-totBodyMass(:,1)),numIt,[])),[nanstd(reshape((totBodyMass(:,2)-totBodyMass(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((totBodyMass(:,2)-totBodyMass(:,1)),numIt,[])))).^0.5)]','r','alpha'); drawnow; set(bl1,'linewidth',2);
% bl1=boundedline([1:numCases], nanmean(reshape((totBodyMass(:,4)-totBodyMass(:,3)),numIt,[])),[nanstd(reshape((totBodyMass(:,4)-totBodyMass(:,3)),numIt,[]))*1.96./((sum(~isnan(reshape((totBodyMass(:,4)-totBodyMass(:,3)),numIt,[])))).^0.5)]','--k','alpha');
% refl=refline(0,0);
% set(refl,'color','k')
% xlim([1 numCases]); xticks([1:numCases]); xticklabels(dispersalLabels);
% xlabel 'dispersal rate'
% ylabel ({'\Delta mean'; 'log_{10}(body size [g])'})
% subplot(4,2,6)
% %violin(reshape(log2(maxBodyMass(:,2)./maxBodyMass(:,1)),numIt,[]),'facecolor','w','facealpha',1,'mc',[],'medc',[]);
% hold on
% bl1=boundedline([1:numCases], nanmean(reshape(consResRatio(:,2)-consResRatio(:,1),numIt,[])),[nanstd(reshape(consResRatio(:,2)-consResRatio(:,1),numIt,[]))*1.96./((sum(~isnan(reshape(consResRatio(:,2)-consResRatio(:,1),numIt,[])))).^0.5)]','r','alpha'); drawnow; set(bl1,'linewidth',2);
% % bl1=boundedline([1:numCases], nanmean(reshape(log2(maxBodyMass(:,2)./maxBodyMass(:,1)),numIt,[])),[nanstd(reshape(log2(maxBodyMass(:,2)./maxBodyMass(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape(log2(maxBodyMass(:,2)./maxBodyMass(:,1)),numIt,[])))).^0.5)]','r','alpha'); drawnow; set(bl1,'linewidth',2);
% % bl1=boundedline([1:numCases], nanmean(reshape(log2(maxBodyMass(:,4)./maxBodyMass(:,3)),numIt,[])),[nanstd(reshape(log2(maxBodyMass(:,4)./maxBodyMass(:,3)),numIt,[]))*1.96./((sum(~isnan(reshape(log2(maxBodyMass(:,4)./maxBodyMass(:,3)),numIt,[])))).^0.5)]','--k','alpha')
% refl=refline(0,0);
% set(refl,'color','k')
% xlim([1 numCases]); xticks([1:numCases]); xticklabels(dispersalLabels);
% xlabel 'dispersal rate'
% % ylabel 'log ratio max body size'
% ylabel '\Delta consumer:resource'


% %plot for one movement rate violin plots of no warming vs. warming in full
% %simulations and in single-species models
% focalMove=2; %focal movement treatment number for the following analyses
% scrsz = get(0,'ScreenSize');
% figs(1)=figure ('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/2 scrsz(4)]);
% subplot(4,2,1)
% violin(log(totBiomass(numIt*(focalMove-1)+1:numIt*focalMove,:)),'facecolor','w','facealpha',1,'mc',[],'medc',[]);
% title ''
% xlabel 'full, full warm, single, single warm'
% ylabel 'ln(biomass)'
% legend off
% subplot(4,2,2)
% hold on
% violin(totProd(numIt*(focalMove-1)+1:numIt*focalMove,:),'facecolor','w','facealpha',1,'mc',[],'medc',[]);
% xlabel 'full, full warm, single, single warm'
% ylabel 'production'
% legend off
% subplot(4,2,3)
% %gscatter([paramIndices;paramIndices],totRich(:),warmingIndices(:),'br','o')
% hold on
% violin(totRich(numIt*(focalMove-1)+1:numIt*focalMove,:),'facecolor','w','facealpha',1,'mc',[],'medc',[]);
% xlabel 'full, full warm, single, single warm'
% ylabel 'richness'
% legend off
% subplot(4,2,4)
% %gscatter([paramIndices;paramIndices],totBeta(:),warmingIndices(:),'br','o')
% hold on
% violin(totBeta(numIt*(focalMove-1)+1:numIt*focalMove,:),'facecolor','w','facealpha',1,'mc',[],'medc',[]);
% xlabel 'full, full warm, single, single warm'
% ylabel 'beta diversity'
% legend off
% subplot(4,2,7)
% %gscatter([paramIndices;paramIndices],totTrophicLevel(:),warmingIndices(:),'br','o')
% hold on
% violin(totTrophicLevel(numIt*(focalMove-1)+1:numIt*focalMove,:),'facecolor','w','facealpha',1,'mc',[],'medc',[]);
% xlabel 'full, full warm'
% ylabel 'mean trophic level'
% legend off
% subplot(4,2,8)
% %gscatter([paramIndices;paramIndices],totTrophicLevel(:),warmingIndices(:),'br','o')
% hold on
% violin(maxTrophicLevel(numIt*(focalMove-1)+1:numIt*focalMove,:),'facecolor','w','facealpha',1,'mc',[],'medc',[]);
% xlabel 'full, full warm'
% ylabel 'max trophic level'
% title ''
% legend off
% subplot(4,2,5)
% %gscatter([paramIndices;paramIndices],totTrophicLevel(:),warmingIndices(:),'br','o')
% hold on
% violin(totBodyMass(numIt*(focalMove-1)+1:numIt*focalMove,:),'facecolor','w','facealpha',1,'mc',[],'medc',[]);
% xlabel 'full, full warm, single, single warm'
% ylabel 'mean body size'
% legend off
% subplot(4,2,6)
% %gscatter([paramIndices;paramIndices],totTrophicLevel(:),warmingIndices(:),'br','o')
% hold on
% violin(maxBodyMass(numIt*(focalMove-1)+1:numIt*focalMove,:),'facecolor','w','facealpha',1,'mc',[],'medc',[]);
% xlabel 'full, full warm, single, single warm'
% ylabel 'max body size'
% title ''
% legend off

% figs(4)=figure ('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/2 3*scrsz(4)/7]);
% subplot(2,2,1)
% hold on
% bl1=boundedline([1:numCases], nanmean(reshape((maxBiomass(:,1)),numIt,[])),[nanstd(reshape((maxBiomass(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((maxBiomass(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2);
% bl1=boundedline([1:numCases], nanmean(reshape((maxBiomass(:,2)),numIt,[])),[nanstd(reshape((maxBiomass(:,2)),numIt,[]))*1.96./((sum(~isnan(reshape((maxBiomass(:,2)),numIt,[])))).^0.5)]','r','alpha'); drawnow; set(bl1,'linewidth',2);
% %bl1=boundedline([1:numCases], nanmean(reshape((maxBiomass(:,3)),numIt,[])),[nanstd(reshape((maxBiomass(:,3)),numIt,[]))*1.96./((sum(~isnan(reshape((maxBiomass(:,3)),numIt,[])))).^0.5)]','--b','alpha','transparency', 0.1);
% %bl1=boundedline([1:numCases], nanmean(reshape((maxBiomass(:,4)),numIt,[])),[nanstd(reshape((maxBiomass(:,4)),numIt,[]))*1.96./((sum(~isnan(reshape((maxBiomass(:,4)),numIt,[])))).^0.5)]','--r','alpha','transparency', 0.1);
% xlabel 'dispersal rate'
% ylabel ({'biomass fraction'; 'from top species'})
% xlim([1 numCases]); xticks([1:numCases]); xticklabels(dispersalLabels);
% title 'outcomes'
% legend off
% subplot(2,2,2)
% hold on
% bl1=boundedline([1:numCases], nanmean(reshape((fractSpeciesProd(:,1)),numIt,[])),[nanstd(reshape((fractSpeciesProd(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((fractSpeciesProd(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2);
% bl1=boundedline([1:numCases], nanmean(reshape((fractSpeciesProd(:,2)),numIt,[])),[nanstd(reshape((fractSpeciesProd(:,2)),numIt,[]))*1.96./((sum(~isnan(reshape((fractSpeciesProd(:,2)),numIt,[])))).^0.5)]','r','alpha'); drawnow; set(bl1,'linewidth',2);
% % bl1=boundedline([1:numCases], nanmean(reshape((maxProd(:,1)),numIt,[])),[nanstd(reshape((maxProd(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((maxProd(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2);
% % bl1=boundedline([1:numCases], nanmean(reshape((maxProd(:,2)),numIt,[])),[nanstd(reshape((maxProd(:,2)),numIt,[]))*1.96./((sum(~isnan(reshape((maxProd(:,2)),numIt,[])))).^0.5)]','r','alpha'); drawnow; set(bl1,'linewidth',2);
% %bl1=boundedline([1:numCases], nanmean(reshape((maxProd(:,3)),numIt,[])),[nanstd(reshape((maxProd(:,3)),numIt,[]))*1.96./((sum(~isnan(reshape((maxProd(:,3)),numIt,[])))).^0.5)]','--b','alpha','transparency', 0.1);
% %bl1=boundedline([1:numCases], nanmean(reshape((maxProd(:,4)),numIt,[])),[nanstd(reshape((maxProd(:,4)),numIt,[]))*1.96./((sum(~isnan(reshape((maxProd(:,4)),numIt,[])))).^0.5)]','--r','alpha','transparency', 0.1);
% xlabel 'dispersal rate'
% ylabel ({'productive fraction'; 'of species'})
% xlim([1 numCases]); xticks([1:numCases]); xticklabels(dispersalLabels);
% title (['warming ' num2str(tempChanges(TempScenario)) '\circC'])
% legend off
% subplot(2,2,3)
% hold on
% bl1=boundedline([1:numCases], nanmean(reshape((maxBiomassBS(:,1)),numIt,[])),[nanstd(reshape((maxBiomassBS(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((maxBiomassBS(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2);
% bl1=boundedline([1:numCases], nanmean(reshape((maxBiomassBS(:,2)),numIt,[])),[nanstd(reshape((maxBiomassBS(:,2)),numIt,[]))*1.96./((sum(~isnan(reshape((maxBiomassBS(:,2)),numIt,[])))).^0.5)]','r','alpha'); drawnow; set(bl1,'linewidth',2);
% %bl1=boundedline([1:numCases], nanmean(reshape((maxBiomassBS(:,3)),numIt,[])),[nanstd(reshape((maxBiomassBS(:,3)),numIt,[]))*1.96./((sum(~isnan(reshape((maxBiomassBS(:,3)),numIt,[])))).^0.5)]','--b','alpha','transparency', 0.1);
% %bl1=boundedline([1:numCases], nanmean(reshape((maxBiomassBS(:,4)),numIt,[])),[nanstd(reshape((maxBiomassBS(:,4)),numIt,[]))*1.96./((sum(~isnan(reshape((maxBiomassBS(:,4)),numIt,[])))).^0.5)]','--r','alpha','transparency', 0.1);
% xlabel 'dispersal rate'
% ylabel ({'most common';'log_{10}(body size [g])'})
% xlim([1 numCases]); xticks([1:numCases]); xticklabels(dispersalLabels);
% title ''
% legend off
% subplot(2,2,4)
% hold on
% bl1=boundedline([1:numCases], nanmean(reshape((maxProdBS(:,1)),numIt,[])),[nanstd(reshape((maxProdBS(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((maxProdBS(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2);
% bl1=boundedline([1:numCases], nanmean(reshape((maxProdBS(:,2)),numIt,[])),[nanstd(reshape((maxProdBS(:,2)),numIt,[]))*1.96./((sum(~isnan(reshape((maxProdBS(:,2)),numIt,[])))).^0.5)]','r','alpha'); drawnow; set(bl1,'linewidth',2);
% %bl1=boundedline([1:numCases], nanmean(reshape((maxProdBS(:,3)),numIt,[])),[nanstd(reshape((maxProdBS(:,3)),numIt,[]))*1.96./((sum(~isnan(reshape((maxProdBS(:,3)),numIt,[])))).^0.5)]','--b','alpha','transparency', 0.1);
% %bl1=boundedline([1:numCases], nanmean(reshape((maxProdBS(:,4)),numIt,[])),[nanstd(reshape((maxProdBS(:,4)),numIt,[]))*1.96./((sum(~isnan(reshape((maxProdBS(:,4)),numIt,[])))).^0.5)]','--r','alpha','transparency', 0.1);
% xlabel 'dispersal rate'
% ylabel ({'most productive';'log_{10}(body size [g])'})
% xlim([1 numCases]); xticks([1:numCases]); xticklabels(dispersalLabels);
% title ''
% legend off
%
% figs(5)=figure ('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/2 3*scrsz(4)/7]);
% subplot(2,2,1)
% hold on
% bl1=boundedline([1:numCases], nanmean(reshape((maxBiomass(:,2)-maxBiomass(:,1)),numIt,[])),[nanstd(reshape((maxBiomass(:,2)-maxBiomass(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((maxBiomass(:,2)-maxBiomass(:,1)),numIt,[])))).^0.5)]','r','alpha'); drawnow; set(bl1,'linewidth',2);
% bl1=boundedline([1:numCases], nanmean(reshape((maxBiomass(:,4)-maxBiomass(:,3)),numIt,[])),[nanstd(reshape((maxBiomass(:,4)-maxBiomass(:,3)),numIt,[]))*1.96./((sum(~isnan(reshape((maxBiomass(:,4)-maxBiomass(:,3)),numIt,[])))).^0.5)]','--k','alpha');
% refl=refline(0,0);
% set(refl,'color','k')
% xlabel 'dispersal rate'
% ylabel ({'\Delta biomass fraction';'from top species'})
% xlim([1 numCases]); xticks([1:numCases]); xticklabels(dispersalLabels);
% title 'warming effects'
% legend off
% subplot(2,2,2)
% hold on
% %bl1=boundedline([1:numCases], nanmean(reshape((maxProd(:,2)-maxProd(:,1)),numIt,[])),[nanstd(reshape((maxProd(:,2)-maxProd(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((maxProd(:,2)-maxProd(:,1)),numIt,[])))).^0.5)]','r','alpha'); drawnow; set(bl1,'linewidth',2);
% %bl1=boundedline([1:numCases], nanmean(reshape((maxProd(:,4)-maxProd(:,3)),numIt,[])),[nanstd(reshape((maxProd(:,4)-maxProd(:,3)),numIt,[]))*1.96./((sum(~isnan(reshape((maxProd(:,4)-maxProd(:,3)),numIt,[])))).^0.5)]','--k','alpha');
% bl1=boundedline([1:numCases], nanmean(reshape((fractSpeciesProd(:,2)-fractSpeciesProd(:,1)),numIt,[])),[nanstd(reshape((fractSpeciesProd(:,2)-fractSpeciesProd(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((fractSpeciesProd(:,2)-fractSpeciesProd(:,1)),numIt,[])))).^0.5)]','r','alpha'); drawnow; set(bl1,'linewidth',2);
% bl1=boundedline([1:numCases], nanmean(reshape((fractSpeciesProd(:,4)-fractSpeciesProd(:,3)),numIt,[])),[nanstd(reshape((fractSpeciesProd(:,4)-fractSpeciesProd(:,3)),numIt,[]))*1.96./((sum(~isnan(reshape((fractSpeciesProd(:,4)-fractSpeciesProd(:,3)),numIt,[])))).^0.5)]','--k','alpha'); drawnow;
% refl=refline(0,0);
% set(refl,'color','k')
% xlabel 'dispersal rate'
% %ylabel '\Delta production from top species'
% ylabel ({'\Delta productive fraction'; 'of species'})
% xlim([1 numCases]); xticks([1:numCases]); xticklabels(dispersalLabels);
% title (['warming ' num2str(tempChanges(TempScenario)) '\circC'])
% legend off
% subplot(2,2,3)
% hold on
% bl1=boundedline([1:numCases], nanmean(reshape((maxBiomassBS(:,2)-maxBiomassBS(:,1)),numIt,[])),[nanstd(reshape((maxBiomassBS(:,2)-maxBiomassBS(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((maxBiomassBS(:,2)-maxBiomassBS(:,1)),numIt,[])))).^0.5)]','r','alpha'); drawnow; set(bl1,'linewidth',2);
% bl1=boundedline([1:numCases], nanmean(reshape((maxBiomassBS(:,4)-maxBiomassBS(:,3)),numIt,[])),[nanstd(reshape((maxBiomassBS(:,4)-maxBiomassBS(:,3)),numIt,[]))*1.96./((sum(~isnan(reshape((maxBiomassBS(:,4)-maxBiomassBS(:,3)),numIt,[])))).^0.5)]','--k','alpha');
% refl=refline(0,0);
% set(refl,'color','k')
% xlabel 'dispersal rate'
% ylabel ({'\Delta most common';'log_{10}(body size [g])'})
% xlim([1 numCases]); xticks([1:numCases]); xticklabels(dispersalLabels);
% title ''
% legend off
% subplot(2,2,4)
% hold on
% bl1=boundedline([1:numCases], nanmean(reshape((maxProdBS(:,2)-maxProdBS(:,1)),numIt,[])),[nanstd(reshape((maxProdBS(:,2)-maxProdBS(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((maxProdBS(:,2)-maxProdBS(:,1)),numIt,[])))).^0.5)]','r','alpha'); drawnow; set(bl1,'linewidth',2);
% bl1=boundedline([1:numCases], nanmean(reshape((maxProdBS(:,4)-maxProdBS(:,3)),numIt,[])),[nanstd(reshape((maxProdBS(:,4)-maxProdBS(:,3)),numIt,[]))*1.96./((sum(~isnan(reshape((maxProdBS(:,4)-maxProdBS(:,3)),numIt,[])))).^0.5)]','--k','alpha');
% refl=refline(0,0);
% set(refl,'color','k')
% xlabel 'dispersal rate'
% ylabel ({'\Delta most productive';'log_{10}(body size [g])'})
% xlim([1 numCases]); xticks([1:numCases]); xticklabels(dispersalLabels);
% title ''
% legend off

figs(6)=figure ('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/1.5 scrsz(4)/1.2]);
subplot(2,2,1)
hold on
bl1=boundedline([1:numCases], nanmean(reshape(rangeShift0(:,1),numIt,[])),[nanstd(reshape(rangeShift0(:,1),numIt,[]))*1.96./((sum(~isnan(reshape(rangeShift0(:,1),numIt,[])))).^0.5)]','r','alpha'); drawnow; set(bl1,'linewidth',2);
bl1=boundedline([1:numCases], nanmean(reshape(rangeShift0(:,2),numIt,[])),[nanstd(reshape(rangeShift0(:,2),numIt,[]))*1.96./((sum(~isnan(reshape(rangeShift0(:,2),numIt,[])))).^0.5)]','--k','alpha');
refl=refline(0,0);
set(refl,'color','k')
%scatter(paramIndices, rangeShift0(:,1),20,'or');
%scatter(paramIndices, rangeShift0(:,2),20,'ok');
scatter(paramIndices+0.5*(rand(size(paramIndices))-0.5), rangeShift0(:,1),'.r');
scatter(paramIndices+0.5*(rand(size(paramIndices))-0.5), rangeShift0(:,2),'.k');
xlabel 'dispersal rate'
ylabel 'means species shift'
title 'no warming'
xlim([0.5 numCases+0.5]); xticks([1:numCases]); xticklabels(dispersalLabels);
legend off

subplot(2,2,2)
hold on
% bl1=boundedline([1:numCases], nanmean(reshape(rangeShift(:,1),numIt,[])),[nanstd(reshape(rangeShift(:,1),numIt,[]))*1.96./((sum(~isnan(reshape(rangeShift(:,1),numIt,[])))).^0.5)]','r','alpha'); drawnow; set(bl1,'linewidth',2);
% bl1=boundedline([1:numCases], nanmean(reshape(rangeShift(:,2),numIt,[])),[nanstd(reshape(rangeShift(:,2),numIt,[]))*1.96./((sum(~isnan(reshape(rangeShift(:,2),numIt,[])))).^0.5)]','--k','alpha');
bl1=boundedline([1:numCases], nanmean(reshape(rangeShift(:,1)-rangeShift0(:,1),numIt,[])),[nanstd(reshape(rangeShift(:,1)-rangeShift0(:,1),numIt,[]))*1.96./((sum(~isnan(reshape(rangeShift(:,1)-rangeShift0(:,1),numIt,[])))).^0.5)]','r','alpha'); drawnow; set(bl1,'linewidth',2);
bl1=boundedline([1:numCases], nanmean(reshape(rangeShift(:,2)-rangeShift0(:,2),numIt,[])),[nanstd(reshape(rangeShift(:,2)-rangeShift0(:,2),numIt,[]))*1.96./((sum(~isnan(reshape(rangeShift(:,2)-rangeShift0(:,2),numIt,[])))).^0.5)]','--k','alpha');
refl=refline(0,0);
set(refl,'color','k')
ref2=refline(0,-3);
set(ref2,'color','r')
% scatter(paramIndices+0.5*(rand(size(paramIndices))-0.5), rangeShift(:,1),'.r');
% scatter(paramIndices+0.5*(rand(size(paramIndices))-0.5), rangeShift(:,2),'.k');
scatter(paramIndices+0.5*(rand(size(paramIndices))-0.5), rangeShift(:,1)-rangeShift0(:,1),'.r');
scatter(paramIndices+0.5*(rand(size(paramIndices))-0.5), rangeShift(:,2)-rangeShift0(:,2),'.k');
xlabel 'dispersal rate'
ylabel 'mean species shift'
xlim([0.5 numCases+0.5]); xticks([1:numCases]); xticklabels(dispersalLabels);
title (['warming-no warming ' num2str(tempChanges(TempScenario)) '\circC'])
legend off

subplot(2,2,3)
hold on
bl1=boundedline([1:numCases], nanmean(reshape(BiomShift0(:,1),numIt,[])),[nanstd(reshape(BiomShift0(:,1),numIt,[]))*1.96./((sum(~isnan(reshape(BiomShift0(:,1),numIt,[])))).^0.5)]','r','alpha'); drawnow; set(bl1,'linewidth',2);
bl1=boundedline([1:numCases], nanmean(reshape(BiomShift0(:,2),numIt,[])),[nanstd(reshape(BiomShift0(:,2),numIt,[]))*1.96./((sum(~isnan(reshape(BiomShift0(:,2),numIt,[])))).^0.5)]','--k','alpha');
refl=refline(0,0);
set(refl,'color','k')
scatter(paramIndices+0.5*(rand(size(paramIndices))-0.5), BiomShift0(:,1),'.r');
scatter(paramIndices+0.5*(rand(size(paramIndices))-0.5), BiomShift0(:,2),'.k');
xlabel 'dispersal rate'
ylabel 'biome shift'
title 'no warming'
xlim([0.5 numCases+0.5]); xticks([1:numCases]); xticklabels(dispersalLabels);
legend off

subplot(2,2,4)
hold on
bl1=boundedline([1:numCases], nanmean(reshape(BiomShift(:,1),numIt,[])),[nanstd(reshape(BiomShift(:,1),numIt,[]))*1.96./((sum(~isnan(reshape(BiomShift(:,1),numIt,[])))).^0.5)]','r','alpha'); drawnow; set(bl1,'linewidth',2);
bl1=boundedline([1:numCases], nanmean(reshape(BiomShift(:,2),numIt,[])),[nanstd(reshape(BiomShift(:,2),numIt,[]))*1.96./((sum(~isnan(reshape(BiomShift(:,2),numIt,[])))).^0.5)]','--k','alpha');
refl=refline(0,0);
set(refl,'color','k')
ref2=refline(0,-3);
set(ref2,'color','r')
scatter(paramIndices+0.5*(rand(size(paramIndices))-0.5), BiomShift(:,1),'.r');
scatter(paramIndices+0.5*(rand(size(paramIndices))-0.5), BiomShift(:,2),'.k');
xlabel 'dispersal rate'
ylabel 'biome shift'
title (['warming ' num2str(tempChanges(TempScenario)) '\circC'])
xlim([0.5 numCases+0.5]); xticks([1:numCases]); xticklabels(dispersalLabels);
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
% pChangeFoodWeb=reshape((totBiomass(:,2)-totBiomass(:,1))./totBiomass(:,1),numIt,[])*100;
% pChangeSingleSp=reshape((totBiomass(:,4)-totBiomass(:,3))./totBiomass(:,3),numIt,[])*100;
pChangeFoodWeb=reshape((totBiomass(:,2)-totBiomass(:,1)),numIt,[]);
pChangeSingleSp=reshape((totBiomass(:,4)-totBiomass(:,3)),numIt,[]);
bl1=boundedline([1:numCases], nanmean(pChangeFoodWeb),[nanstd(pChangeFoodWeb)*1.96./((sum(~isnan(pChangeFoodWeb))).^0.5)]','--r','alpha'); drawnow; set(bl1,'linewidth',2);
bl1=boundedline([1:numCases], nanmean(pChangeSingleSp),[nanstd(pChangeSingleSp)*1.96./((sum(~isnan(pChangeSingleSp))).^0.5)]','--k','alpha');
refl=refline(0,0);
set(refl,'color','r')
xlim([1 numCases]); xticks([1:numCases]); xticklabels(dispersalLabels);
%title 'outcomes and warming effects'
title ''
xlabel 'dispersal rate'
ylabel '\Delta'
legend off
subplot(4,2,2)
yyaxis left
hold on
bl1=boundedline([1:numCases], nanmean(reshape((totProd(:,1)),numIt,[])),[nanstd(reshape((totProd(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((totProd(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2);
ylabel 'production [gm^-{3}/day]'
yyaxis right
hold on
% pChangeFoodWeb=reshape((totProd(:,2)-totProd(:,1))./totProd(:,1),numIt,[])*100;
% pChangeSingleSp=reshape((totProd(:,4)-totProd(:,3))./totProd(:,3),numIt,[])*100;
pChangeFoodWeb=reshape((totProd(:,2)-totProd(:,1)),numIt,[]);
pChangeSingleSp=reshape((totProd(:,4)-totProd(:,3)),numIt,[]);
bl1=boundedline([1:numCases], nanmean(pChangeFoodWeb),[nanstd(pChangeFoodWeb)*1.96./((sum(~isnan(pChangeFoodWeb))).^0.5)]','--r','alpha'); drawnow; set(bl1,'linewidth',2);
bl1=boundedline([1:numCases], nanmean(pChangeSingleSp),[nanstd(pChangeSingleSp)*1.96./((sum(~isnan(pChangeSingleSp))).^0.5)]','--k','alpha');
refl=refline(0,0);
set(refl,'color','r')
xlim([1 numCases]); xticks([1:numCases]); xticklabels(dispersalLabels);
title (['warming ' num2str(tempChanges(TempScenario)) '\circC'])
xlabel 'dispersal rate'
ylabel '\Delta'
legend off
subplot(4,2,3)
yyaxis left
hold on
bl1=boundedline([1:numCases], nanmean(reshape((totRich(:,1)),numIt,[])),[nanstd(reshape((totRich(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((totRich(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2);
ylabel 'richness'
yyaxis right
hold on
% pChangeFoodWeb=reshape((totRich(:,2)-totRich(:,1))./totRich(:,1),numIt,[])*100;
% pChangeSingleSp=reshape((totRich(:,4)-totRich(:,3))./totRich(:,3),numIt,[])*100;
pChangeFoodWeb=reshape((totRich(:,2)-totRich(:,1)),numIt,[]);
pChangeSingleSp=reshape((totRich(:,4)-totRich(:,3)),numIt,[]);
bl1=boundedline([1:numCases], nanmean(pChangeFoodWeb),[nanstd(pChangeFoodWeb)*1.96./((sum(~isnan(pChangeFoodWeb))).^0.5)]','--r','alpha'); drawnow; set(bl1,'linewidth',2);
bl1=boundedline([1:numCases], nanmean(pChangeSingleSp),[nanstd(pChangeSingleSp)*1.96./((sum(~isnan(pChangeSingleSp))).^0.5)]','--k','alpha');
refl=refline(0,0);
set(refl,'color','r')
xlim([1 numCases]); xticks([1:numCases]); xticklabels(dispersalLabels);
xlabel 'dispersal rate'
ylabel '\Delta'
title ''
legend off
subplot(4,2,4)
yyaxis left
hold on
bl1=boundedline([1:numCases], nanmean(reshape((totBeta(:,1)),numIt,[])),[nanstd(reshape((totBeta(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((totBeta(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2);
ylabel 'beta diversity'
yyaxis right
% pChangeFoodWeb=reshape((totBeta(:,2)-totBeta(:,1))./totBeta(:,1),numIt,[])*100;
% pChangeSingleSp=reshape((totBeta(:,4)-totBeta(:,3))./totBeta(:,3),numIt,[])*100;
pChangeFoodWeb=reshape((totBeta(:,2)-totBeta(:,1)),numIt,[]);
pChangeSingleSp=reshape((totBeta(:,4)-totBeta(:,3)),numIt,[]);
bl1=boundedline([1:numCases], nanmean(pChangeFoodWeb),[nanstd(pChangeFoodWeb)*1.96./((sum(~isnan(pChangeFoodWeb))).^0.5)]','--r','alpha'); drawnow; set(bl1,'linewidth',2);
bl1=boundedline([1:numCases], nanmean(pChangeSingleSp),[nanstd(pChangeSingleSp)*1.96./((sum(~isnan(pChangeSingleSp))).^0.5)]','--k','alpha');
refl=refline(0,0);
set(refl,'color','r')
xlim([1 numCases]); xticks([1:numCases]); xticklabels(dispersalLabels);
xlabel 'dispersal rate'
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
%pChangeFoodWeb=reshape((totTrophicLevel(:,2)-totTrophicLevel(:,1))./totTrophicLevel(:,1),numIt,[])*100;
pChangeFoodWeb=reshape((totTrophicLevel(:,2)-totTrophicLevel(:,1)),numIt,[]);
bl1=boundedline([1:numCases], nanmean(pChangeFoodWeb),[nanstd(pChangeFoodWeb)*1.96./((sum(~isnan(pChangeFoodWeb))).^0.5)]','--r','alpha'); drawnow; set(bl1,'linewidth',2);
refl=refline(0,0);
set(refl,'color','r')
xlim([1 numCases]); xticks([1:numCases]); xticklabels(dispersalLabels);
xlabel 'dispersal rate'
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
%pChangeFoodWeb=reshape((maxTrophicLevel(:,2)-maxTrophicLevel(:,1))./maxTrophicLevel(:,1),numIt,[])*100;
pChangeFoodWeb=reshape((maxTrophicLevel(:,2)-maxTrophicLevel(:,1)),numIt,[]);
bl1=boundedline([1:numCases], nanmean(pChangeFoodWeb),[nanstd(pChangeFoodWeb)*1.96./((sum(~isnan(pChangeFoodWeb))).^0.5)]','--r','alpha'); drawnow; set(bl1,'linewidth',2);
refl=refline(0,0);
set(refl,'color','r')
xlim([1 numCases]); xticks([1:numCases]); xticklabels(dispersalLabels);
xlabel 'dispersal rate'
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
% pChangeFoodWeb=reshape((totBodyMass(:,2)-totBodyMass(:,1))./totBodyMass(:,1),numIt,[])*100;
% pChangeSingleSp=reshape((totBodyMass(:,4)-totBodyMass(:,3))./totBodyMass(:,3),numIt,[])*100;
pChangeFoodWeb=reshape((totBodyMass(:,2)-totBodyMass(:,1)),numIt,[]);
pChangeSingleSp=reshape((totBodyMass(:,4)-totBodyMass(:,3)),numIt,[]);
bl1=boundedline([1:numCases], nanmean(pChangeFoodWeb),[nanstd(pChangeFoodWeb)*1.96./((sum(~isnan(pChangeFoodWeb))).^0.5)]','--r','alpha'); drawnow; set(bl1,'linewidth',2);
bl1=boundedline([1:numCases], nanmean(pChangeSingleSp),[nanstd(pChangeSingleSp)*1.96./((sum(~isnan(pChangeSingleSp))).^0.5)]','--k','alpha');
refl=refline(0,0);
set(refl,'color','r')
xlim([1 numCases]); xticks([1:numCases]); xticklabels(dispersalLabels);
xlabel 'dispersal rate'
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
%pChangeFoodWeb=reshape((consResRatio(:,2)-consResRatio(:,1))./consResRatio(:,1),numIt,[])*100;
pChangeFoodWeb=reshape((consResRatio(:,2)-consResRatio(:,1)),numIt,[]);
bl1=boundedline([1:numCases], nanmean(pChangeFoodWeb),[nanstd(pChangeFoodWeb)*1.96./((sum(~isnan(pChangeFoodWeb))).^0.5)]','--r','alpha'); drawnow; set(bl1,'linewidth',2);
refl=refline(0,0);
set(refl,'color','r')
xlim([1 numCases]); xticks([1:numCases]); xticklabels(dispersalLabels);
xlabel 'dispersal rate'
ylabel '\Delta'
title ''
legend off



compfig=figure ('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/1.75 3*scrsz(4)/7]);
set(compfig,'defaultAxesColorOrder',[[0 0 1]; [1 0 0]]);
subplot(2,2,1)
yyaxis left
bl1=boundedline([1:numCases], nanmean(reshape((maxBiomass(:,1)),numIt,[])),[nanstd(reshape((maxBiomass(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((maxBiomass(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2);
ylabel ({'biomass fraction'; 'from top species'})
yyaxis right
hold on
% pChangeFoodWeb=reshape((maxBiomass(:,2)-maxBiomass(:,1))./maxBiomass(:,1),numIt,[])*100;
% pChangeSingleSp=reshape((maxBiomass(:,4)-maxBiomass(:,3))./maxBiomass(:,3),numIt,[])*100;
pChangeFoodWeb=reshape((maxBiomass(:,2)-maxBiomass(:,1)),numIt,[]);
pChangeSingleSp=reshape((maxBiomass(:,4)-maxBiomass(:,3)),numIt,[]);
bl1=boundedline([1:numCases], nanmean(pChangeFoodWeb),[nanstd(pChangeFoodWeb)*1.96./((sum(~isnan(pChangeFoodWeb))).^0.5)]','--r','alpha'); drawnow; set(bl1,'linewidth',2);
bl1=boundedline([1:numCases], nanmean(pChangeSingleSp),[nanstd(pChangeSingleSp)*1.96./((sum(~isnan(pChangeSingleSp))).^0.5)]','--k','alpha');
refl=refline(0,0);
set(refl,'color','r')
ylabel '\Delta'
xlabel 'dispersal rate'
xlim([1 numCases]); xticks([1:numCases]); xticklabels(dispersalLabels);
%title 'outcomes'
title ''
legend off
subplot(2,2,2)
yyaxis left
bl1=boundedline([1:numCases], nanmean(reshape((fractSpeciesProd(:,1)),numIt,[])),[nanstd(reshape((fractSpeciesProd(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((fractSpeciesProd(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2);
ylabel ({'productive fraction'; 'of species'})
yyaxis right
hold on
% pChangeFoodWeb=reshape((fractSpeciesProd(:,2)-fractSpeciesProd(:,1))./fractSpeciesProd(:,1),numIt,[])*100;
% pChangeSingleSp=reshape((fractSpeciesProd(:,4)-fractSpeciesProd(:,3))./fractSpeciesProd(:,3),numIt,[])*100;
pChangeFoodWeb=reshape((fractSpeciesProd(:,2)-fractSpeciesProd(:,1)),numIt,[]);
pChangeSingleSp=reshape((fractSpeciesProd(:,4)-fractSpeciesProd(:,3)),numIt,[]);
bl1=boundedline([1:numCases], nanmean(pChangeFoodWeb),[nanstd(pChangeFoodWeb)*1.96./((sum(~isnan(pChangeFoodWeb))).^0.5)]','--r','alpha'); drawnow; set(bl1,'linewidth',2);
bl1=boundedline([1:numCases], nanmean(pChangeSingleSp),[nanstd(pChangeSingleSp)*1.96./((sum(~isnan(pChangeSingleSp))).^0.5)]','--k','alpha');
refl=refline(0,0);
set(refl,'color','r')
xlabel 'dispersal rate'
ylabel '\Delta'
xlim([1 numCases]); xticks([1:numCases]); xticklabels(dispersalLabels);
title (['warming ' num2str(tempChanges(TempScenario)) '\circC'])
legend off
subplot(2,2,3)
yyaxis left
bl1=boundedline([1:numCases], nanmean(reshape((maxBiomassBS(:,1)),numIt,[])),[nanstd(reshape((maxBiomassBS(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((maxBiomassBS(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2);
ylabel ({'most common';'log_{10}(body size [g])'})
yyaxis right
hold on
% pChangeFoodWeb=reshape((maxBiomassBS(:,2)-maxBiomassBS(:,1))./maxBiomassBS(:,1),numIt,[])*100;
% pChangeSingleSp=reshape((maxBiomassBS(:,4)-maxBiomassBS(:,3))./maxBiomassBS(:,3),numIt,[])*100;
pChangeFoodWeb=reshape((maxBiomassBS(:,2)-maxBiomassBS(:,1)),numIt,[]);
pChangeSingleSp=reshape((maxBiomassBS(:,4)-maxBiomassBS(:,3)),numIt,[]);
bl1=boundedline([1:numCases], nanmean(pChangeFoodWeb),[nanstd(pChangeFoodWeb)*1.96./((sum(~isnan(pChangeFoodWeb))).^0.5)]','--r','alpha'); drawnow; set(bl1,'linewidth',2);
bl1=boundedline([1:numCases], nanmean(pChangeSingleSp),[nanstd(pChangeSingleSp)*1.96./((sum(~isnan(pChangeSingleSp))).^0.5)]','--k','alpha');
refl=refline(0,0);
set(refl,'color','r')
xlabel 'dispersal rate'
ylabel '\Delta'
xlim([1 numCases]); xticks([1:numCases]); xticklabels(dispersalLabels);
title ''
legend off
subplot(2,2,4)
yyaxis left
bl1=boundedline([1:numCases], nanmean(reshape((maxProdBS(:,1)),numIt,[])),[nanstd(reshape((maxProdBS(:,1)),numIt,[]))*1.96./((sum(~isnan(reshape((maxProdBS(:,1)),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(bl1,'linewidth',2);
ylabel ({'most productive';'log_{10}(body size [g])'})
yyaxis right
hold on
% pChangeFoodWeb=reshape((maxProdBS(:,2)-maxProdBS(:,1))./maxProdBS(:,1),numIt,[])*100;
% pChangeSingleSp=reshape((maxProdBS(:,4)-maxProdBS(:,3))./maxProdBS(:,3),numIt,[])*100;
pChangeFoodWeb=reshape((maxProdBS(:,2)-maxProdBS(:,1)),numIt,[]);
pChangeSingleSp=reshape((maxProdBS(:,4)-maxProdBS(:,3)),numIt,[]);
bl1=boundedline([1:numCases], nanmean(pChangeFoodWeb),[nanstd(pChangeFoodWeb)*1.96./((sum(~isnan(pChangeFoodWeb))).^0.5)]','--r','alpha'); drawnow; set(bl1,'linewidth',2);
bl1=boundedline([1:numCases], nanmean(pChangeSingleSp),[nanstd(pChangeSingleSp)*1.96./((sum(~isnan(pChangeSingleSp))).^0.5)]','--k','alpha');
refl=refline(0,0);
set(refl,'color','k')
xlabel 'dispersal rate'
ylabel '\Delta'
xlim([1 numCases]); xticks([1:numCases]); xticklabels(dispersalLabels);
title ''
legend off


shiftfig=figure ('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/2 scrsz(4)/1.5]); %3.5
set(shiftfig,'defaultAxesColorOrder',[[0 0 1]; [1 0 0]]);
subplot(2,2,2) %median shift
BS_move=reshape(allBody_Biomass(1,:,:),numIt*size(allBody_Biomass,2),[]); %reshape individual species body sizes by movement rates
Shift_move=reshape(All_Shifts(2,:,:),numIt*size(All_Shifts,2),[]); %reshape individual species shifts with warming by movement rates
meanBSShiftCoeff=[];
seBSShiftCoeff=[];
meanBS26Shift=[];
seBS26Shift=[];
for i=1:length(moveRates)
    lm=fitlm(BS_move(:,i),Shift_move(:,i));
    %meanBSShiftCoeff(:,i)=lm.Coefficients.Estimate;
    %seBSShiftCoeff(:,i)=lm.Coefficients.SE;
    [shift2, ci2] = predict(lm,2);
    [shift6, ci6] = predict(lm,6);
    meanBS26Shift(:,i)=[shift2;shift6];
    seBS26Shift(:,i)=[ci2(1);ci6(1)];
end
hold on
%yyaxis left
CM=colormap(jet(128)); % set colormap
%bl1=boundedline([1:numCases], meanBSShiftCoeff(2,:),seBSShiftCoeff(2,:)*1.96,'-b','alpha'); drawnow; set(bl1,'linewidth',2);
bl2=boundedline([1:numCases], meanBS26Shift(1,:),[meanBS26Shift(1,:)-seBS26Shift(1,:)],'--','Cmap',CM(ceil(2*128/6),:),'alpha'); drawnow; set(bl2,'linewidth',1);
bl6=boundedline([1:numCases], meanBS26Shift(2,:),[meanBS26Shift(2,:)-seBS26Shift(2,:)],'--','Cmap',CM(end,:),'alpha'); drawnow; set(bl6,'linewidth',1);
if ~exist('CentroidShift1') %if this is first set of analysis (specialists)
    plotColor='b';
else
    plotColor='r';
end
blmean=boundedline([1:numCases], nanmean(reshape(rangeShift(:,1)-rangeShift0(:,1),numIt,[])),[nanstd(reshape(rangeShift(:,1)-rangeShift0(:,1),numIt,[]))*1.96./((sum(~isnan(reshape(rangeShift(:,1)-rangeShift0(:,1),numIt,[])))).^0.5)]',plotColor,'alpha'); drawnow; set(blmean,'linewidth',4);
blssproj=boundedline([1:numCases], nanmean(reshape(rangeShift(:,2)-rangeShift0(:,2),numIt,[])),[nanstd(reshape(rangeShift(:,2)-rangeShift0(:,2),numIt,[]))*1.96./((sum(~isnan(reshape(rangeShift(:,2)-rangeShift0(:,2),numIt,[])))).^0.5)]','--k','alpha'); set(blssproj,'linewidth',2);
refl=refline(0,0);
set(refl,'color','k')
ylabel 'median shift'
ylim([-3.2 0])
%yyaxis right
refl=refline(0,-3);
set(refl,'color',plotColor)
%ylabel 'log_{10}(body size):shift intercept'
xlim([1 numCases-1]); xticks([1:numCases-1]); xticklabels(dispersalLabels);
xlabel 'dispersal rate'

subplot(2,2,1) %centroid shift
CentroidShift_move=reshape(All_CentroidShifts(2,:,:),numIt*size(All_CentroidShifts,2),[]); %reshape individual species shifts with warming by movement rates
meanBSShiftCoeff=[];
seBSShiftCoeff=[];
meanBS26Shift=[];
seBS26Shift=[];
for i=1:length(moveRates)
    lm=fitlm(BS_move(:,i),CentroidShift_move(:,i));
    %meanBSShiftCoeff(:,i)=lm.Coefficients.Estimate;
    %seBSShiftCoeff(:,i)=lm.Coefficients.SE;
    [shift2, ci2] = predict(lm,2);
    [shift6, ci6] = predict(lm,6);
    meanBS26Shift(:,i)=[shift2;shift6];
    seBS26Shift(:,i)=[ci2(1);ci6(1)];
end
hold on
%yyaxis left
CM=colormap(jet(128)); % set colormap
%bl1=boundedline([1:numCases], meanBSShiftCoeff(2,:),seBSShiftCoeff(2,:)*1.96,'-b','alpha'); drawnow; set(bl1,'linewidth',2);
bl2=boundedline([1:numCases], meanBS26Shift(1,:),[meanBS26Shift(1,:)-seBS26Shift(1,:)],'--','Cmap',CM(ceil(2*128/6),:),'alpha'); drawnow; set(bl2,'linewidth',1);
bl6=boundedline([1:numCases], meanBS26Shift(2,:),[meanBS26Shift(2,:)-seBS26Shift(2,:)],'--','Cmap',CM(end,:),'alpha'); drawnow; set(bl6,'linewidth',1);
if ~exist('CentroidShift1') %if this is first set of analysis (specialists)
    plotColor='b';
else
    plotColor='r';
end
blmean=boundedline([1:numCases], nanmean(reshape(CentroidShift(:,1)-CentroidShift0(:,1),numIt,[])),[nanstd(reshape(CentroidShift(:,1)-CentroidShift0(:,1),numIt,[]))*1.96./((sum(~isnan(reshape(CentroidShift(:,1)-CentroidShift0(:,1),numIt,[])))).^0.5)]',plotColor,'alpha'); drawnow; set(blmean,'linewidth',4);
blssproj=boundedline([1:numCases], nanmean(reshape(CentroidShift(:,2)-CentroidShift0(:,2),numIt,[])),[nanstd(reshape(CentroidShift(:,2)-CentroidShift0(:,2),numIt,[]))*1.96./((sum(~isnan(reshape(CentroidShift(:,2)-CentroidShift0(:,2),numIt,[])))).^0.5)]','--k','alpha'); set(blssproj,'linewidth',2);
refl=refline(0,0);
set(refl,'color','k')
ylabel 'centroid shift'
ylim([-3.2 0])
%yyaxis right
refl=refline(0,-3);
set(refl,'color',plotColor)
%ylabel 'log_{10}(body size):shift intercept'
xlim([1 numCases-1]); xticks([1:numCases-1]); xticklabels(dispersalLabels);
xlabel 'dispersal rate'

subplot(2,2,3) %leading edge shift
LeadingShift_move=reshape(All_LeadingShifts(2,:,:),numIt*size(All_LeadingShifts,2),[]); %reshape individual species shifts with warming by movement rates
meanBSShiftCoeff=[];
seBSShiftCoeff=[];
meanBS26Shift=[];
seBS26Shift=[];
for i=1:length(moveRates)
    lm=fitlm(BS_move(:,i),LeadingShift_move(:,i));
    %meanBSShiftCoeff(:,i)=lm.Coefficients.Estimate;
    %seBSShiftCoeff(:,i)=lm.Coefficients.SE;
    [shift2, ci2] = predict(lm,2);
    [shift6, ci6] = predict(lm,6);
    meanBS26Shift(:,i)=[shift2;shift6];
    seBS26Shift(:,i)=[ci2(1);ci6(1)];
end
hold on
%yyaxis left
CM=colormap(jet(128)); % set colormap
%bl1=boundedline([1:numCases], meanBSShiftCoeff(2,:),seBSShiftCoeff(2,:)*1.96,'-b','alpha'); drawnow; set(bl1,'linewidth',2);
bl2=boundedline([1:numCases], meanBS26Shift(1,:),[meanBS26Shift(1,:)-seBS26Shift(1,:)],'--','Cmap',CM(ceil(2*128/6),:),'alpha'); drawnow; set(bl2,'linewidth',1);
bl6=boundedline([1:numCases], meanBS26Shift(2,:),[meanBS26Shift(2,:)-seBS26Shift(2,:)],'--','Cmap',CM(end,:),'alpha'); drawnow; set(bl6,'linewidth',1);
if ~exist('CentroidShift1') %if this is first set of analysis (specialists)
    plotColor='b';
else
    plotColor='r';
end
blmean=boundedline([1:numCases], nanmean(reshape(LeadingShift(:,1)-LeadingShift0(:,1),numIt,[])),[nanstd(reshape(LeadingShift(:,1)-LeadingShift0(:,1),numIt,[]))*1.96./((sum(~isnan(reshape(LeadingShift(:,1)-LeadingShift0(:,1),numIt,[])))).^0.5)]',plotColor,'alpha'); drawnow; set(blmean,'linewidth',4);
blssproj=boundedline([1:numCases], nanmean(reshape(LeadingShift(:,2)-LeadingShift0(:,2),numIt,[])),[nanstd(reshape(LeadingShift(:,2)-LeadingShift0(:,2),numIt,[]))*1.96./((sum(~isnan(reshape(LeadingShift(:,2)-LeadingShift0(:,2),numIt,[])))).^0.5)]','--k','alpha'); set(blssproj,'linewidth',2);
refl=refline(0,0);
set(refl,'color','k')
ylabel 'leading edge shift'
ylim([-3.2 0])
%yyaxis right
refl=refline(0,-3);
set(refl,'color',plotColor)
%ylabel 'log_{10}(body size):shift intercept'
xlim([1 numCases-1]); xticks([1:numCases-1]); xticklabels(dispersalLabels);
xlabel 'dispersal rate'

subplot(2,2,4) %trailing edge shift
TrailingShift_move=reshape(All_TrailingShifts(2,:,:),numIt*size(All_TrailingShifts,2),[]); %reshape individual species shifts with warming by movement rates
meanBSShiftCoeff=[];
seBSShiftCoeff=[];
meanBS26Shift=[];
seBS26Shift=[];
for i=1:length(moveRates)
    lm=fitlm(BS_move(:,i),TrailingShift_move(:,i));
    %meanBSShiftCoeff(:,i)=lm.Coefficients.Estimate;
    %seBSShiftCoeff(:,i)=lm.Coefficients.SE;
    [shift2, ci2] = predict(lm,2);
    [shift6, ci6] = predict(lm,6);
    meanBS26Shift(:,i)=[shift2;shift6];
    seBS26Shift(:,i)=[ci2(1);ci6(1)];
end
hold on
%yyaxis left
CM=colormap(jet(128)); % set colormap
%bl1=boundedline([1:numCases], meanBSShiftCoeff(2,:),seBSShiftCoeff(2,:)*1.96,'-b','alpha'); drawnow; set(bl1,'linewidth',2);
bl2=boundedline([1:numCases], meanBS26Shift(1,:),[meanBS26Shift(1,:)-seBS26Shift(1,:)],'--','Cmap',CM(ceil(2*128/6),:),'alpha'); drawnow; set(bl2,'linewidth',1);
bl6=boundedline([1:numCases], meanBS26Shift(2,:),[meanBS26Shift(2,:)-seBS26Shift(2,:)],'--','Cmap',CM(end,:),'alpha'); drawnow; set(bl6,'linewidth',1);
if ~exist('CentroidShift1') %if this is first set of analysis (specialists)
    plotColor='b';
else
    plotColor='r';
end
blmean=boundedline([1:numCases], nanmean(reshape(TrailingShift(:,1)-TrailingShift0(:,1),numIt,[])),[nanstd(reshape(TrailingShift(:,1)-TrailingShift0(:,1),numIt,[]))*1.96./((sum(~isnan(reshape(TrailingShift(:,1)-TrailingShift0(:,1),numIt,[])))).^0.5)]',plotColor,'alpha'); drawnow; set(blmean,'linewidth',4);
blssproj=boundedline([1:numCases], nanmean(reshape(TrailingShift(:,2)-TrailingShift0(:,2),numIt,[])),[nanstd(reshape(TrailingShift(:,2)-TrailingShift0(:,2),numIt,[]))*1.96./((sum(~isnan(reshape(TrailingShift(:,2)-TrailingShift0(:,2),numIt,[])))).^0.5)]','--k','alpha'); set(blssproj,'linewidth',2);
refl=refline(0,0);
set(refl,'color','k')
ylabel 'trailing edge shift'
ylim([-3.2 0])
%yyaxis right
refl=refline(0,-3);
set(refl,'color',plotColor)
%ylabel 'log_{10}(body size):shift intercept'
xlim([1 numCases-1]); xticks([1:numCases-1]); xticklabels(dispersalLabels);
xlabel 'dispersal rate'

if ~exist('CentroidShift1') %record percent change in food web species shift compared to single-species projection
    %percShift1=reshape((rangeShift(:,1)-rangeShift0(:,1)-(rangeShift(:,2)-rangeShift0(:,2)))./(rangeShift(:,2)-rangeShift0(:,2)),numIt,[])*100;
    CentroidShift1=CentroidShift;
    CentroidShift01=CentroidShift0;
    rangeShift1=rangeShift;
    rangeShift01=rangeShift0;
    TrailingShift1=TrailingShift;
    TrailingShift01=TrailingShift0;
    LeadingShift1=LeadingShift;
    LeadingShift01=LeadingShift0;
    %percShift1=reshape((rangeShift(:,1)-rangeShift0(:,1)-(rangeShift(:,2)-rangeShift0(:,2)))./(rangeShift(:,2)-rangeShift0(:,2)),numIt,[])*100;
    %percShift1((rangeShift(:,1)-rangeShift0(:,1)-(rangeShift(:,2)-rangeShift0(:,2)))==0)=0; % (% difference in shift) is 0 when the absolute shift is 0, even when the expected shift (demoninator in the line above) is 0.
    
else %plot both sets of changes
    shiftfig2=figure ('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/2 scrsz(4)/1.5]); %3.5
    set(shiftfig2,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);
    subplot(2,2,2) %median shift
    %yyaxis left
    hold on
    blmean2=boundedline([1:numCases], nanmean(reshape(rangeShift(:,1)-rangeShift0(:,1),numIt,[])),[nanstd(reshape(rangeShift(:,1)-rangeShift0(:,1),numIt,[]))*1.96./((sum(~isnan(reshape(rangeShift(:,1)-rangeShift0(:,1),numIt,[])))).^0.5)]','r','alpha'); drawnow; set(blmean2,'linewidth',4);
    blmean1=boundedline([1:numCases], nanmean(reshape(rangeShift1(:,1)-rangeShift01(:,1),numIt,[])),[nanstd(reshape(rangeShift1(:,1)-rangeShift01(:,1),numIt,[]))*1.96./((sum(~isnan(reshape(rangeShift1(:,1)-rangeShift01(:,1),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(blmean1,'linewidth',4);
    %blssproj=boundedline([1:numCases], nanmean(reshape(rangeShift(:,2)-rangeShift0(:,2),numIt,[])),[nanstd(reshape(rangeShift(:,2)-rangeShift0(:,2),numIt,[]))*1.96./((sum(~isnan(reshape(rangeShift(:,2)-rangeShift0(:,2),numIt,[])))).^0.5)]','--k','alpha'); set(blssproj,'linewidth',2);
    refl=refline(0,0);
    set(refl,'color','k')
    ylabel 'median shift'
    ylim([-3.2 0])
%     yyaxis right
%     hold on
%     percShift1=reshape((rangeShift1(:,1)-rangeShift01(:,1)-(rangeShift1(:,2)-rangeShift01(:,2)))./(rangeShift1(:,2)-rangeShift01(:,2)),numIt,[])*100;
%     percShift1((rangeShift1(:,1)-rangeShift01(:,1)-(rangeShift1(:,2)-rangeShift01(:,2)))==0)=0; % (% difference in shift) is 0 when the absolute shift is 0, even when the expected shift (demoninator in the line above) is 0.
%     percShift2=reshape((rangeShift(:,1)-rangeShift0(:,1)-(rangeShift(:,2)-rangeShift0(:,2)))./(rangeShift(:,2)-rangeShift0(:,2)),numIt,[])*100;
%     percShift2((rangeShift(:,1)-rangeShift0(:,1)-(rangeShift(:,2)-rangeShift0(:,2)))==0)=0; % (% difference in shift) is 0 when the absolute shift is 0, even when the expected shift (demoninator in the line above) is 0.
%     blpercShift1=boundedline([1:numCases], nanmean(percShift1),[nanstd(percShift1)*1.96./((sum(~isnan(percShift1))).^0.5)]','--b','alpha'); drawnow; set(blpercShift1,'linewidth',2); %if running specialists first (plotted blue), then generalists (red)
%     blpercShift2=boundedline([1:numCases], nanmean(percShift2),[nanstd(percShift2)*1.96./((sum(~isnan(percShift2))).^0.5)]','--r','alpha'); drawnow; set(blpercShift2,'linewidth',2);
%     ylabel '% median shift'
    refl=refline(0,-3);
    set(refl,'color',plotColor)
    %ylabel 'log_{10}(body size):shift intercept'
    xlim([1 numCases-1]); xticks([1:numCases-1]); xticklabels(dispersalLabels);
    xlabel 'dispersal rate'
    
    subplot(2,2,1) %centroid shift
    %yyaxis left
    hold on
    blmean2=boundedline([1:numCases], nanmean(reshape(CentroidShift(:,1)-CentroidShift0(:,1),numIt,[])),[nanstd(reshape(CentroidShift(:,1)-CentroidShift0(:,1),numIt,[]))*1.96./((sum(~isnan(reshape(CentroidShift(:,1)-CentroidShift0(:,1),numIt,[])))).^0.5)]','r','alpha'); drawnow; set(blmean2,'linewidth',4);
    blmean1=boundedline([1:numCases], nanmean(reshape(CentroidShift1(:,1)-CentroidShift01(:,1),numIt,[])),[nanstd(reshape(CentroidShift1(:,1)-CentroidShift01(:,1),numIt,[]))*1.96./((sum(~isnan(reshape(CentroidShift1(:,1)-CentroidShift01(:,1),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(blmean1,'linewidth',4);
    %blssproj=boundedline([1:numCases], nanmean(reshape(CentroidShift(:,2)-CentroidShift0(:,2),numIt,[])),[nanstd(reshape(CentroidShift(:,2)-CentroidShift0(:,2),numIt,[]))*1.96./((sum(~isnan(reshape(CentroidShift(:,2)-CentroidShift0(:,2),numIt,[])))).^0.5)]','--k','alpha'); set(blssproj,'linewidth',2);
    refl=refline(0,0);
    set(refl,'color','k')
    ylabel 'centroid shift'
    ylim([-3.2 0])
%     yyaxis right
%     hold on
%     percShift1=reshape((CentroidShift1(:,1)-CentroidShift01(:,1)-(CentroidShift1(:,2)-CentroidShift01(:,2)))./(CentroidShift1(:,2)-CentroidShift01(:,2)),numIt,[])*100;
%     percShift1((CentroidShift1(:,1)-CentroidShift01(:,1)-(CentroidShift1(:,2)-CentroidShift01(:,2)))==0)=0; % (% difference in shift) is 0 when the absolute shift is 0, even when the expected shift (demoninator in the line above) is 0.
%     percShift2=reshape((CentroidShift(:,1)-CentroidShift0(:,1)-(CentroidShift(:,2)-CentroidShift0(:,2)))./(CentroidShift(:,2)-CentroidShift0(:,2)),numIt,[])*100;
%     percShift2((CentroidShift(:,1)-CentroidShift0(:,1)-(CentroidShift(:,2)-CentroidShift0(:,2)))==0)=0; % (% difference in shift) is 0 when the absolute shift is 0, even when the expected shift (demoninator in the line above) is 0.
%     blpercShift1=boundedline([1:numCases], nanmean(percShift1),[nanstd(percShift1)*1.96./((sum(~isnan(percShift1))).^0.5)]','--b','alpha'); drawnow; set(blpercShift1,'linewidth',2); %if running specialists first (plotted blue), then generalists (red)
%     blpercShift2=boundedline([1:numCases], nanmean(percShift2),[nanstd(percShift2)*1.96./((sum(~isnan(percShift2))).^0.5)]','--r','alpha'); drawnow; set(blpercShift2,'linewidth',2);
%     ylabel '% centroid shift'
    refl=refline(0,-3);
    set(refl,'color',plotColor)
    %ylabel 'log_{10}(body size):shift intercept'
    xlim([1 numCases-1]); xticks([1:numCases-1]); xticklabels(dispersalLabels);
    xlabel 'dispersal rate'
    
    subplot(2,2,3) %leading edge shift
    %yyaxis left
    hold on
    blmean2=boundedline([1:numCases], nanmean(reshape(LeadingShift(:,1)-LeadingShift0(:,1),numIt,[])),[nanstd(reshape(LeadingShift(:,1)-LeadingShift0(:,1),numIt,[]))*1.96./((sum(~isnan(reshape(LeadingShift(:,1)-LeadingShift0(:,1),numIt,[])))).^0.5)]','r','alpha'); drawnow; set(blmean2,'linewidth',4);
    blmean1=boundedline([1:numCases], nanmean(reshape(LeadingShift1(:,1)-LeadingShift01(:,1),numIt,[])),[nanstd(reshape(LeadingShift1(:,1)-LeadingShift01(:,1),numIt,[]))*1.96./((sum(~isnan(reshape(LeadingShift1(:,1)-LeadingShift01(:,1),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(blmean1,'linewidth',4);
    %blssproj=boundedline([1:numCases], nanmean(reshape(LeadingShift(:,2)-LeadingShift0(:,2),numIt,[])),[nanstd(reshape(LeadingShift(:,2)-LeadingShift0(:,2),numIt,[]))*1.96./((sum(~isnan(reshape(LeadingShift(:,2)-LeadingShift0(:,2),numIt,[])))).^0.5)]','--k','alpha'); set(blssproj,'linewidth',2);
    refl=refline(0,0);
    set(refl,'color','k')
    ylabel 'leading edge shift'
    ylim([-3.2 0])
%     yyaxis right
%     hold on
%     percShift1=reshape((LeadingShift1(:,1)-LeadingShift01(:,1)-(LeadingShift1(:,2)-LeadingShift01(:,2)))./(LeadingShift1(:,2)-LeadingShift01(:,2)),numIt,[])*100;
%     percShift1((LeadingShift1(:,1)-LeadingShift01(:,1)-(LeadingShift1(:,2)-LeadingShift01(:,2)))==0)=0; % (% difference in shift) is 0 when the absolute shift is 0, even when the expected shift (demoninator in the line above) is 0.
%     percShift2=reshape((LeadingShift(:,1)-LeadingShift0(:,1)-(LeadingShift(:,2)-LeadingShift0(:,2)))./(LeadingShift(:,2)-LeadingShift0(:,2)),numIt,[])*100;
%     percShift2((LeadingShift(:,1)-LeadingShift0(:,1)-(LeadingShift(:,2)-LeadingShift0(:,2)))==0)=0; % (% difference in shift) is 0 when the absolute shift is 0, even when the expected shift (demoninator in the line above) is 0.
%     blpercShift1=boundedline([1:numCases], nanmean(percShift1),[nanstd(percShift1)*1.96./((sum(~isnan(percShift1))).^0.5)]','--b','alpha'); drawnow; set(blpercShift1,'linewidth',2); %if running specialists first (plotted blue), then generalists (red)
%     blpercShift2=boundedline([1:numCases], nanmean(percShift2),[nanstd(percShift2)*1.96./((sum(~isnan(percShift2))).^0.5)]','--r','alpha'); drawnow; set(blpercShift2,'linewidth',2);
%     ylabel '% leading edge shift'
    refl=refline(0,-3);
    set(refl,'color',plotColor)
    %ylabel 'log_{10}(body size):shift intercept'
    xlim([1 numCases-1]); xticks([1:numCases-1]); xticklabels(dispersalLabels);
    xlabel 'dispersal rate'
    
    subplot(2,2,4) %trailing edge shift
    %yyaxis left
    blmean2=boundedline([1:numCases], nanmean(reshape(TrailingShift(:,1)-TrailingShift0(:,1),numIt,[])),[nanstd(reshape(TrailingShift(:,1)-TrailingShift0(:,1),numIt,[]))*1.96./((sum(~isnan(reshape(TrailingShift(:,1)-TrailingShift0(:,1),numIt,[])))).^0.5)]','r','alpha'); drawnow; set(blmean2,'linewidth',4);
    blmean1=boundedline([1:numCases], nanmean(reshape(TrailingShift1(:,1)-TrailingShift01(:,1),numIt,[])),[nanstd(reshape(TrailingShift1(:,1)-TrailingShift01(:,1),numIt,[]))*1.96./((sum(~isnan(reshape(TrailingShift1(:,1)-TrailingShift01(:,1),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(blmean1,'linewidth',4);
    %blssproj=boundedline([1:numCases], nanmean(reshape(TrailingShift(:,2)-TrailingShift0(:,2),numIt,[])),[nanstd(reshape(TrailingShift(:,2)-TrailingShift0(:,2),numIt,[]))*1.96./((sum(~isnan(reshape(TrailingShift(:,2)-TrailingShift0(:,2),numIt,[])))).^0.5)]','--k','alpha'); set(blssproj,'linewidth',2);
    refl=refline(0,0);
    set(refl,'color','k')
    ylabel 'trailing edge shift'
    ylim([-3.2 0])
%     yyaxis right
%     hold on
%     percShift1=reshape((TrailingShift1(:,1)-TrailingShift01(:,1)-(TrailingShift1(:,2)-TrailingShift01(:,2)))./(TrailingShift1(:,2)-TrailingShift01(:,2)),numIt,[])*100;
%     percShift1((TrailingShift1(:,1)-TrailingShift01(:,1)-(TrailingShift1(:,2)-TrailingShift01(:,2)))==0)=0; % (% difference in shift) is 0 when the absolute shift is 0, even when the expected shift (demoninator in the line above) is 0.
%     percShift2=reshape((TrailingShift(:,1)-TrailingShift0(:,1)-(TrailingShift(:,2)-TrailingShift0(:,2)))./(TrailingShift(:,2)-TrailingShift0(:,2)),numIt,[])*100;
%     percShift2((TrailingShift(:,1)-TrailingShift0(:,1)-(TrailingShift(:,2)-TrailingShift0(:,2)))==0)=0; % (% difference in shift) is 0 when the absolute shift is 0, even when the expected shift (demoninator in the line above) is 0.
%     blpercShift1=boundedline([1:numCases], nanmean(percShift1),[nanstd(percShift1)*1.96./((sum(~isnan(percShift1))).^0.5)]','--b','alpha'); drawnow; set(blpercShift1,'linewidth',2); %if running specialists first (plotted blue), then generalists (red)
%     blpercShift2=boundedline([1:numCases], nanmean(percShift2),[nanstd(percShift2)*1.96./((sum(~isnan(percShift2))).^0.5)]','--r','alpha'); drawnow; set(blpercShift2,'linewidth',2);
%     ylabel '% trailing edge shift'
    refl=refline(0,-3);
    set(refl,'color',plotColor)
    %ylabel 'log_{10}(body size):shift intercept'
    xlim([1 numCases-1]); xticks([1:numCases-1]); xticklabels(dispersalLabels);
    xlabel 'dispersal rate'
    %
    %     shiftfig2=figure ('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/1.5 scrsz(4)/2.5]); %3.5
    %     set(shiftfig2,'defaultAxesColorOrder',[[0 0 1]; [1 0 0]]);
    %     subplot(1,2,1)
    %     hold on
    %     percShift2=reshape((rangeShift(:,1)-rangeShift0(:,1)-(rangeShift(:,2)-rangeShift0(:,2)))./(rangeShift(:,2)-rangeShift0(:,2)),numIt,[])*100;
    %     precShift2((rangeShift(:,1)-rangeShift0(:,1)-(rangeShift(:,2)-rangeShift0(:,2)))==0)=0; % (% difference in shift) is 0 when the absolute shift is 0, even when the expected shift (demoninator in the line above) is 0.
    %     blpercShift1=boundedline([1:numCases], nanmean(percShift1),[nanstd(percShift1)*1.96./((sum(~isnan(percShift1))).^0.5)]','b','alpha'); drawnow; set(blpercShift1,'linewidth',4); %if running specialists first (plotted blue), then generalists (red)
    %     blpercShift2=boundedline([1:numCases], nanmean(percShift2),[nanstd(percShift2)*1.96./((sum(~isnan(percShift2))).^0.5)]','r','alpha'); drawnow; set(blpercShift2,'linewidth',4);
    %     refl=refline(0,0);
    %     set(refl,'color','k')
    %     xlim([2 numCases-1]); xticks([2:numCases-1]); xticklabels(moveRates(2:end));
    %     xlabel 'dispersal rate'
    %     ylabel '% difference in median shift'
    
    shiftfig3=figure ('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/3 scrsz(4)/3]); %3.5
    set(shiftfig3,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);
    hold on
    %plot(nanmean(reshape(rangeShift(:,1)-rangeShift0(:,1),numIt,[])),'.-r','linewidth',2);
    %plot(nanmean(reshape(rangeShift1(:,1)-rangeShift01(:,1),numIt,[])),'.-b','linewidth',2);
    refl=refline(0,0);
    set(refl,'color','k')
    blLVCentroid2=boundedline([1:numCases], nanmean(reshape(CentroidShift(:,2)-CentroidShift0(:,2),numIt,[])),[nanstd(reshape(CentroidShift(:,2)-CentroidShift0(:,2),numIt,[]))*1.96./((sum(~isnan(reshape(CentroidShift(:,2)-CentroidShift0(:,2),numIt,[])))).^0.5)]','k','alpha'); drawnow; set(blLVCentroid2,'linewidth',2);
    blLVCentroid1=boundedline([1:numCases], nanmean(reshape(CentroidShift1(:,2)-CentroidShift01(:,2),numIt,[])),[nanstd(reshape(CentroidShift1(:,2)-CentroidShift01(:,2),numIt,[]))*1.96./((sum(~isnan(reshape(CentroidShift1(:,2)-CentroidShift01(:,2),numIt,[])))).^0.5)]','k','alpha'); drawnow; set(blLVCentroid1,'linewidth',2);
    blCentroid2=boundedline([1:numCases], nanmean(reshape(CentroidShift(:,1)-CentroidShift0(:,1),numIt,[])),[nanstd(reshape(CentroidShift(:,1)-CentroidShift0(:,1),numIt,[]))*1.96./((sum(~isnan(reshape(CentroidShift(:,1)-CentroidShift0(:,1),numIt,[])))).^0.5)]','r','alpha'); drawnow; set(blCentroid2,'linewidth',4);
    blCentroid1=boundedline([1:numCases], nanmean(reshape(CentroidShift1(:,1)-CentroidShift01(:,1),numIt,[])),[nanstd(reshape(CentroidShift1(:,1)-CentroidShift01(:,1),numIt,[]))*1.96./((sum(~isnan(reshape(CentroidShift1(:,1)-CentroidShift01(:,1),numIt,[])))).^0.5)]','b','alpha'); drawnow; set(blCentroid1,'linewidth',4);
    plot([1:numCases], nanmean(reshape(LeadingShift(:,1)-LeadingShift0(:,1),numIt,[])),'--r','linewidth',2);
    plot([1:numCases], nanmean(reshape(LeadingShift1(:,1)-LeadingShift01(:,1),numIt,[])),'--b','linewidth',2);
    plot([1:numCases], nanmean(reshape(TrailingShift(:,1)-TrailingShift0(:,1),numIt,[])),'--r','linewidth',2);
    plot([1:numCases], nanmean(reshape(TrailingShift1(:,1)-TrailingShift01(:,1),numIt,[])),'--b','linewidth',2);
%     lowerDist2=nanmean(reshape(CentroidShift(:,1)-CentroidShift0(:,1),numIt,[]))-nanmean(reshape(TrailingShift(:,1)-TrailingShift0(:,1),numIt,[]));
%     upperDist2=nanmean(reshape(LeadingShift(:,1)-LeadingShift0(:,1),numIt,[]))-nanmean(reshape(CentroidShift(:,1)-CentroidShift0(:,1),numIt,[]));
%     lowerDist1=nanmean(reshape(CentroidShift1(:,1)-CentroidShift01(:,1),numIt,[]))-nanmean(reshape(TrailingShift(:,1)-TrailingShift0(:,1),numIt,[]));
%     upperDist1=nanmean(reshape(TrailingShift1(:,1)-TrailingShift01(:,1),numIt,[]))-nanmean(reshape(CentroidShift1(:,1)-CentroidShift01(:,1),numIt,[]));
%     blleadtrail2=boundedline([1:numCases], nanmean(reshape(CentroidShift(:,1)-CentroidShift0(:,1),numIt,[])),lowerDist2,upperDist2,'--r','alpha'); drawnow; set(blleadtrail2,'linewidth',2);
%     blleadtrail1=boundedline([1:numCases], nanmean(reshape(CentroidShift1(:,1)-CentroidShift01(:,1),numIt,[])),lowerDist1,upperDist1,'--b','alpha'); drawnow; set(blleadtrail1,'linewidth',2);
    refl=refline(0,-3);
    set(refl,'color',plotColor)
    ylabel 'range shift'
    ylim([-3.2 0])
    xlim([1 numCases-1]); xticks([1:numCases-1]); xticklabels(dispersalLabels);
    xlabel 'dispersal rate'
end

nfig=figure ('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/4.2 scrsz(4)/4]);
foodwebNoWarmingSurvive=sum(reshape(totRich(:,1)>0,numIt,[]));
foodwebWarmingSurvive=sum(reshape(totRich(:,2)>0,numIt,[]));
singlespeciesSurvive=sum(reshape(totRich(:,3)>0,numIt,[]));
hold on
bar([1:numCases],singlespeciesSurvive,'w')
bar([1:numCases],foodwebNoWarmingSurvive,'b')
bar([1:numCases],foodwebWarmingSurvive,'r')
xlim([0 numCases+1]); xticks([1:numCases]); xticklabels(dispersalLabels);
xlabel 'dispersal rate'
ylabel 'non-extinct food webs'
