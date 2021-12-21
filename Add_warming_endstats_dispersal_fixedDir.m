
%Edward Wong Oct 22, 17
%collect data from Make_warming_MultPtsstats_Parallel0 simulations

clear

TimeData=string(datetime);

%choose temperature change scenario to display: (position of tempChange =[-2 2 4 6])
TempScenario=1; %position is without 0 change position [+2 +4 +6] (3 for base case or 1 for variants)
recordYrs=20; %use these number of years at the end of time series to compute mean statistics
defaultNumIt=10; %cap maximum number of iterations per case to analyse

%specify simulation data folder in current directory
Path='Type III narrow thermal evlp'; %this is an example

%%
%use these lines for specialist food webs:
%AntiCases={'basalSize0.01_meanD-Inf_pInedible0_fIII';'basalSize0.01_meanD0_pInedible0_fIII';'basalSize0.01_meanD3_pInedible0_fIII';'basalSize0.01_meanD6_pInedible0_fIII';'basalSize0.01_meanD7_pInedible0_fIII';'basalSize0.01_meanD8_pInedible0_fIII';'basalSize0.01_meanD9_pInedible0_fIII';'basalSize0.01_meanD10_pInedible0_fIII'};
%Cases={'basalSize0.01_meanD-Inf_pInedible0.5_fIII';'basalSize0.01_meanD0_pInedible0.5_fIII';'basalSize0.01_meanD3_pInedible0.5_fIII';'basalSize0.01_meanD6_pInedible0.5_fIII';'basalSize0.01_meanD7_pInedible0.5_fIII';'basalSize0.01_meanD8_pInedible0.5_fIII';'basalSize0.01_meanD9_pInedible0.5_fIII';'basalSize0.01_meanD10_pInedible0.5_fIII'};

%use these lines for generalist food webs:
Cases={'basalSize0.01_meanD-Inf_pInedible0_fIII';'basalSize0.01_meanD0_pInedible0_fIII';'basalSize0.01_meanD3_pInedible0_fIII';'basalSize0.01_meanD6_pInedible0_fIII';'basalSize0.01_meanD7_pInedible0_fIII';'basalSize0.01_meanD8_pInedible0_fIII';'basalSize0.01_meanD9_pInedible0_fIII';'basalSize0.01_meanD10_pInedible0_fIII'};
AntiCases={'basalSize0.01_meanD-Inf_pInedible0.5_fIII';'basalSize0.01_meanD0_pInedible0.5_fIII';'basalSize0.01_meanD3_pInedible0.5_fIII';'basalSize0.01_meanD6_pInedible0.5_fIII';'basalSize0.01_meanD7_pInedible0.5_fIII';'basalSize0.01_meanD8_pInedible0.5_fIII';'basalSize0.01_meanD9_pInedible0.5_fIII';'basalSize0.01_meanD10_pInedible0.5_fIII'};
%%

numCases=length(Cases);
moveRates=[-Inf 0 3 6 7 8]; %diffusion rates
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
novelSpatialAssemblage=[]; %fraction of current spatial occurances that is new
lostSpatialAssemblage=[]; %fraction of original spatial occurances lost
novelSpatialAssemblage2=[]; %fraction of current spatial occurances that is new
lostSpatialAssemblage2=[]; %fraction of original spatial occurances lost
novelSpatialAssemblage6=[]; %fraction of current spatial occurances that is new
lostSpatialAssemblage6=[]; %fraction of original spatial occurances lost
novelCoexistence=[]; %fraction of current coexisting pairs that is new
lostCoexistence=[]; %fraction of originial coexisting pairs lost
novelCoexistence22=[]; %fraction of current coexisting pairs that is new
lostCoexistence22=[]; %fraction of originial coexisting pairs lost
novelCoexistence26=[]; %fraction of current coexisting pairs that is new
lostCoexistence26=[]; %fraction of originial coexisting pairs lost
novelCoexistence66=[]; %fraction of current coexisting pairs that is new
lostCoexistence66=[]; %fraction of originial coexisting pairs lost

%by space
varBiomass=[];
varProd=[];
varTrophicLevel=[];
%by body mass
allBody_Biomass=[]; %record all non-extinct species. col1: body mass, col2: opt temp, col3: biomass for no warming cases, col4: biomass for warming cases.
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
Centroid2Shift0=[]; %centroid shift with no warming
Centroid2Shift=[]; %centroid shift with warming
Centroid6Shift0=[]; %centroid shift with no warming
Centroid6Shift=[]; %centroid shift with warming
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
RangeExpansion0=[]; %average range expansion under no warming
RangeExpansion=[]; %average range expansion under warming
Trailing2Shift0=[]; %average trailing edge shift under no warming
Trailing2Shift=[]; %average trailing edge shift under warming
Leading2Shift0=[]; %average trailing edge shift under no warming
Leading2Shift=[]; %average trailing edge shift under warming
Range2Expansion0=[]; %average range expansion under no warming
Range2Expansion=[]; %average range expansion under warmingTrailingShift0=[]; %average trailing edge shift under no warming
Trailing6Shift0=[]; %average trailing edge shift under no warming
Trailing6Shift=[]; %average trailing edge shift under warming
Leading6Shift0=[]; %average trailing edge shift under no warming
Leading6Shift=[]; %average trailing edge shift under warming
Range6Expansion0=[]; %average range expansion under no warming
Range6Expansion=[]; %average range expansion under warming
RangeExpansionPerc0=[]; %average range expansion % under no warming
RangeExpansionPerc=[]; %average range expansion % under warming
Range2ExpansionPerc0=[]; %average range expansion % under no warming
Range2ExpansionPerc=[]; %average range expansion % under warming
Range6ExpansionPerc0=[]; %average range expansion % under no warming
Range6ExpansionPerc=[]; %average range expansion % under warming


%species-specific variables
All_BodySize=[]; %record all surviving species body sizes in all scenarios
All_OptTemp=[]; %optimal temperature
All_InitTemp=[]; %initial patch temperature
All_Move=[]; %movement rate
All_TempChange=[]; %temperature change over 200 years
All_ShiftLag=[]; %lag in species shift
All_Shifts=[]; %all median location cls
All_CentroidShifts=[]; %all mean (centroid) shifts
All_TrailingShifts=[];
All_LeadingShifts=[];
All_RangeExpansions=[];
All_RangeSizes=[]; %all initial range sizes
All_Centroids=[]; %all initial centroids

All_a_est=[];
All_r_est=[];

% 'Select files containing simulation results in first folder'
%
% Files={};
% Paths={};
%
% reply1 = 'a';
%
% while (reply1 == 'a')
%[FileName,PathName]=uigetfile('.mat','Simulation data','MultiSelect','on')
%newFiles=cellstr(FileName);
%newPaths=cellstr(PathName);
%Files=[Files, newFiles];
AllFiles=struct2cell(dir(Path));
Files=AllFiles(1,:);
% for f=1:length(newFiles)
%     Paths=[Paths, newPaths];
% end
% reply1 = input('Finished? (enter to continue, a to add files, q to quit): ','s');
% if isempty(reply1)
%     reply1 = 'continue';
% end
%end

numFiles=length(Files);

%if (reply1 ~='q') %continue
for CaseNumber=1:numCases %first, check for minimum number of iterations for each case
    iteration=0;
    Positions=contains(Files,Cases{CaseNumber}); %find positions of .mat files that belong to movement rate treatment of CaseNumber
    AntiPositions=contains(Files,AntiCases{CaseNumber});
    for filePos=1:numFiles
        if(Positions(filePos)==1 && AntiPositions(filePos)~=1)
            iteration=iteration+1;
            caseIts(CaseNumber)=iteration; %update number of iterations for each case
        end
    end
end
numIt=min([caseIts defaultNumIt]); %minimum number of iterations for all cases

numRun=1;
for CaseNumber=1:numCases
    iteration=0;
    Positions=contains(Files,Cases{CaseNumber}); %find positions of .mat files that belong to landscape type CaseNumber
    AntiPositions=contains(Files,AntiCases{CaseNumber});
    for filePos=1:numFiles
        if(Positions(filePos)==1 && AntiPositions(filePos)~=1 && iteration<numIt)
            iteration=iteration+1;
            load([Path '/' Files{filePos}]);
            
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
            
            tempChanges=P.dT*73000; %temperature changes corresponding to temperature scenario
            TempChange=tempChanges(TempScenario);
            
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
            
            %Redo single-species LV hindcast fits by fixing c and estimating r and a for each species:
            display(['running single species estimate for run ' num2str(numRun) ' of ' num2str(numIt*numCases)])
            numRun=numRun+1;
            [r5,a5,c5,z5,K5,flag5,raR25,r_T5,K_T5,K_T_ratio5,r_T_ratio5]=estSingleSpeciesModelmsy(Btrans,dBtrans,gainBtrans,P,[0.5 0.5]); %fit growth model to all patches at once using estimated best quantiles
            All_r_est=[All_r_est; r5];
            All_a_est=[All_a_est; a5];
            %                                 BPquantiles=[BPquantiles;fitCodeEsts];
            %                                 BPpercDiff=[BPpercDiff;[percDiffB,percDiffP]];
            %                                 disp(['case ' num2str(CaseNumber) ' iteration ' num2str(iteration) ': quantiles ' num2str(fitCodeEsts,2) ', %diff ' num2str([percDiffB,percDiffP],2)])
            %
            %rerun single-species projection:
            TimePts=[1:365:200*365+1]; %record every year
            numPts=length(TimePts); %number of time points
            gainBLV5=zeros(P.nx,P.n);
            dBLV5=zeros(P.nx,P.n);
            gainBLV5w=zeros(P.nx,P.n);
            dBLV5w=zeros(P.nx,P.n);
            BLV5=BLV1_yrs(:,:,1); %no warming case under estimated single-species dynamics (start at identical point as other simulations)
            BLV5w=BLV5; %6C warming case under estimated single-species dynamics (start at identical point as other simulations)
            for t = 1:TimePts(end) %run for 200 years (with daily time steps)
                BLV5(BLV5<eps) = 0;
                BLV5w(BLV5w<eps)= 0;
                T1      = P.T;% + t.*P.dT; %<<< add this when time is right
                T1w      = P.T + (t-1)*P.dT;
                BLV5_RK1=sub_move(BLV5,P); %single species model no temp change Runge-Kutta step 1
                [gainBLV5_RK1 dBLV5_RK1] = sub_demogLV(BLV5_RK1,T1,r5,a5,c5,z5,P.Ea,P.k,P.s.mi,P.Spd);
                BLV5_RK2        = sub_move(BLV5+0.5*dBLV5_RK1,P); % move Runge-Kutta step 2
                [gainBLV5_RK2 dBLV5_RK2] = sub_demogLV(BLV5_RK2,T1,r5,a5,c5,z5,P.Ea,P.k,P.s.mi,P.Spd);
                dBLV5=(dBLV5_RK1+dBLV5_RK2)/2; %2nd order RK integration
                gainBLV5=(gainBLV5_RK1+gainBLV5_RK2)/2; %2nd order RK integration
                
                BLV5w_RK1=sub_move(BLV5w,P); %single species model temp change Runge-Kutta step 1
                [gainBLV5w_RK1 dBLV5w_RK1] = sub_demogLV(BLV5w_RK1,T1w,r5,a5,c5,z5,P.Ea,P.k,P.s.mi,P.Spd);
                BLV5w_RK2        = sub_move(BLV5w+0.5*dBLV5w_RK1,P); % move Runge-Kutta step 2
                [gainBLV5w_RK2 dBLV5w_RK2] = sub_demogLV(BLV5w_RK2,T1w,r5,a5,c5,z5,P.Ea,P.k,P.s.mi,P.Spd);
                dBLV5w=(dBLV5w_RK1+dBLV5w_RK2)/2; %2nd order RK integration
                gainBLV5w=(gainBLV5w_RK1+gainBLV5w_RK2)/2; %2nd order RK integration
                tpos=find(t==TimePts);
                if ~isempty(tpos)
                    BLV5_yrs(:,:,tpos)=BLV5;
                    gainBLV5_yrs(:,:,tpos)=gainBLV5;
                    BLV5w_yrs(:,:,tpos)=BLV5w;
                    gainBLV5w_yrs(:,:,tpos)=gainBLV5w;
                end
                BLV5 = BLV5 + dBLV5;
                BLV5w = BLV5w + dBLV5w;
            end
            BLV1=nanmean(BLV5_yrs(:,:,end-recordYrs+1:end),3);
            gainBLV1=nanmean(gainBLV5_yrs(:,:,end-recordYrs+1:end),3);
            BLV1w=nanmean(BLV5w_yrs(:,:,end-recordYrs+1:end),3);
            gainBLV1w=nanmean(gainBLV5w_yrs(:,:,end-recordYrs+1:end),3);
            
            
            %---------------Change assignments here for different single-species projections---------
            %             BLV1=nanmean(BLV1_yrs(:,:,end-recordYrs+1:end),3); %BLV, BLV1-6
            %             gainBLV1=nanmean(gainBLV1_yrs(:,:,end-recordYrs+1:end),3);
            %             BLV1w=nanmean(BLV1w_yrs(:,:,end-recordYrs+1:end),3);
            %             gainBLV1w=nanmean(gainBLV1w_yrs(:,:,end-recordYrs+1:end),3); %gainBLV1w=nanmean(gainBLV4w_yrs(:,:,end-recordYrs+1:end,TempScenario),3);
            %----------------------------------------------------------------------------------------
            %end
            %    end
            BtransEnd=nanmean(Btrans(:,:,end-recordYrs+1:end),3); %biomasses at the end of transcient period
            BtransEnd(BtransEnd<=eps)=nan;
            numPatches=size(B,1);
            %-------
            % %             %%use if eliminating species with centroids in the coldest patches:
            %             meanInitCentroid=nansum(BtransEnd.*[1:numPatches]')./nansum(BtransEnd); %get initial centroids
            %             BtransEnd(:,meanInitCentroid<=3)=NaN; %set initial biomass of these species to NaN
            %%use if eliminating species with leading edges in 3 coldest patches and trailing edges in 3 hottest patches:
            leadingQ=0.025; %lower quantile indicates location of leading edge (towards cold region)
            trailingQ=0.975; %upper quantile indicates location of trailing edge (towards hot region)
            minSpeciesBiomass=leadingQ*nansum(BtransEnd); %determine leading biomass density
            maxSpeciesBiomass=trailingQ*nansum(BtransEnd); %determine trailing biomass density
            tempInitTrailing=sum(cumsum(BtransEnd,'omitnan')<maxSpeciesBiomass)+1;
            tempInitLeading=sum(cumsum(BtransEnd,'omitnan')<minSpeciesBiomass)+1;
            BtransEnd(:,tempInitLeading<=3)=NaN; %set initial biomass of these species to NaN
            B(:,tempInitLeading<=3)=NaN; %set initial biomass of these species to NaN
            Bw(:,tempInitLeading<=3)=NaN; %set initial biomass of these species to NaN
            BLV1(:,tempInitLeading<=3)=NaN; %set initial biomass of these species to NaN
            BLV1w(:,tempInitLeading<=3)=NaN; %set initial biomass of these species to NaN
            gainB(:,tempInitLeading<=3)=NaN;
            gainBw(:,tempInitLeading<=3)=NaN;
            gainBLV1(:,tempInitLeading<=3)=NaN;
            gainBLV1w(:,tempInitLeading<=3)=NaN;
            
            %-------
            
            Bnan=B;
            Bnan(Bnan<=eps)=nan;
            BLV1nan=BLV1;
            BLV1nan(BLV1nan<=eps)=nan;
            Bwnan=Bw;
            Bwnan(Bwnan<=eps)=nan;
            BLV1wnan=BLV1w;
            BLV1wnan(BLV1wnan<=eps)=nan;
            
            %variables for no warming/warming (col1, col2 unless otherwise commented) cases:
            %global (for heterotrophs)
            totBiomass=[totBiomass; [mean(nansum(B,2)) mean(nansum(Bw,2)) mean(nansum(BLV1,2)) mean(nansum(BLV1w,2))]]; %total biomass (as mean biomass/m^3)
            allBody_Biomass=cat(3,allBody_Biomass,[P.S(2:end);P.z;nansum(B);nansum(Bw);nansum(BLV1);nansum(BLV1w)]);
            totProd=[totProd; [mean(nansum(gainB,2)) mean(nansum(gainBw,2)) mean(nansum(gainBLV1,2)) mean(nansum(gainBLV1w,2))]]; %total production (as mean production/(day m^3))
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
            BtransExist=BtransEnd>eps; %matrix of species existence (columns) at different locations (rows) at the end of transcient period
            Bexist=Bnan>eps; %matrix of food web species existence (columns) at different locations (rows) at the end
            Bwexist=Bwnan>eps;
            BLV1exist=BLV1nan>eps;
            BLV1wexist=BLV1wnan>eps;
            BLocalDiff=Bexist-BtransExist; %new (+1) or lost (-1) local species
            BwLocalDiff=Bwexist-BtransExist;
            BLV1LocalDiff=BLV1exist-BtransExist;
            BLV1wLocalDiff=BLV1wexist-BtransExist;
            novelSpatialAssemblage=[novelSpatialAssemblage; [sum(BLocalDiff(:)==1)/sum(Bexist(:)) sum(BwLocalDiff(:)==1)/sum(Bwexist(:)) sum(BLV1LocalDiff(:)==1)/sum(BLV1exist(:)) sum(BLV1wLocalDiff(:)==1)/sum(BLV1wexist(:))]]; % portion of current species that were not there initially
            lostSpatialAssemblage=[lostSpatialAssemblage; [sum(BLocalDiff(:)==-1)/sum(BtransExist(:)) sum(BwLocalDiff(:)==-1)/sum(BtransExist(:)) sum(BLV1LocalDiff(:)==-1)/sum(BtransExist(:)) sum(BLV1wLocalDiff(:)==-1)/sum(BtransExist(:))]]; %portion of initially present species that are no longer there
            Size2_3_Pos=find(P.S>2&P.S<=3)-1;
            Size5_6_Pos=find(P.S>5&P.S<=6)-1;
            novelSpatialAssemblage2=[novelSpatialAssemblage2; [sum(sum(BLocalDiff(:,Size2_3_Pos)==1))/sum(sum(Bexist(:,Size2_3_Pos))) sum(sum(BwLocalDiff(:,Size2_3_Pos)==1))/sum(sum(Bexist(:,Size2_3_Pos))) sum(sum(BLV1LocalDiff(:,Size2_3_Pos)==1))/sum(sum(BLV1exist(:,Size2_3_Pos))) sum(sum(BLV1wLocalDiff(:,Size2_3_Pos)==1))/sum(sum(BLV1exist(:,Size2_3_Pos)))]]; % portion of current species that were not there initially
            novelSpatialAssemblage6=[novelSpatialAssemblage6; [sum(sum(BLocalDiff(:,Size5_6_Pos)==1))/sum(sum(Bexist(:,Size5_6_Pos))) sum(sum(BwLocalDiff(:,Size5_6_Pos)==1))/sum(sum(Bexist(:,Size5_6_Pos))) sum(sum(BLV1LocalDiff(:,Size5_6_Pos)==1))/sum(sum(BLV1exist(:,Size5_6_Pos))) sum(sum(BLV1wLocalDiff(:,Size5_6_Pos)==1))/sum(sum(BLV1exist(:,Size5_6_Pos)))]]; % portion of current species that were not there initially
            lostSpatialAssemblage2=[lostSpatialAssemblage2; [sum(sum(BLocalDiff(:,Size2_3_Pos)==-1))/sum(sum(BtransExist(:,Size2_3_Pos))) sum(sum(BwLocalDiff(:,Size2_3_Pos)==-1))/sum(sum(BtransExist(:,Size2_3_Pos))) sum(sum(BLV1LocalDiff(:,Size2_3_Pos)==-1))/sum(sum(BtransExist(:,Size2_3_Pos))) sum(sum(BLV1wLocalDiff(:,Size2_3_Pos)==-1))/sum(sum(BtransExist(:,Size2_3_Pos)))]]; %portion of initially present species that are no longer there
            lostSpatialAssemblage6=[lostSpatialAssemblage6; [sum(sum(BLocalDiff(:,Size5_6_Pos)==-1))/sum(sum(BtransExist(:,Size5_6_Pos))) sum(sum(BwLocalDiff(:,Size5_6_Pos)==-1))/sum(sum(BtransExist(:,Size5_6_Pos))) sum(sum(BLV1LocalDiff(:,Size5_6_Pos)==-1))/sum(sum(BtransExist(:,Size5_6_Pos))) sum(sum(BLV1wLocalDiff(:,Size5_6_Pos)==-1))/sum(sum(BtransExist(:,Size5_6_Pos)))]]; %portion of initially present species that are no longer there
            
            
            for sp1=1:length(B)
                for sp2=1:length(B)
                    CoexistInit(sp1,sp2)=sum(BtransEnd(:,sp1).*BtransEnd(:,sp2)>0)>0; %initial coexistence matrix
                    CoexistB(sp1,sp2)=sum(B(:,sp1).*B(:,sp2)>0)>0; %final coexistence matrix
                    CoexistBw(sp1,sp2)=sum(Bw(:,sp1).*Bw(:,sp2)>0)>0; %final coexistence matrix
                    CoexistBLV1(sp1,sp2)=sum(BLV1(:,sp1).*BLV1(:,sp2)>0)>0; %final coexistence matrix
                    CoexistBLV1w(sp1,sp2)=sum(BLV1w(:,sp1).*BLV1w(:,sp2)>0)>0; %final coexistence matrix
                    if (~isempty(find(Size2_3_Pos==sp1)))
                        if (~isempty(find(Size2_3_Pos==sp2))) %within size class coexistence
                            Coexist22Init(sp1,sp2)=sum(BtransEnd(:,sp1).*BtransEnd(:,sp2)>0)>0; %initial coexistence matrix
                            Coexist22B(sp1,sp2)=sum(B(:,sp1).*B(:,sp2)>0)>0; %final coexistence matrix
                            Coexist22Bw(sp1,sp2)=sum(Bw(:,sp1).*Bw(:,sp2)>0)>0; %final coexistence matrix
                            Coexist22BLV1(sp1,sp2)=sum(BLV1(:,sp1).*BLV1(:,sp2)>0)>0; %final coexistence matrix
                            Coexist22BLV1w(sp1,sp2)=sum(BLV1w(:,sp1).*BLV1w(:,sp2)>0)>0; %final coexistence matrix
                        elseif (~isempty(find(Size5_6_Pos==sp2))) %between size class coexistence
                            Coexist26Init(sp1,sp2)=sum(BtransEnd(:,sp1).*BtransEnd(:,sp2)>0)>0; %initial coexistence matrix
                            Coexist26B(sp1,sp2)=sum(B(:,sp1).*B(:,sp2)>0)>0; %final coexistence matrix
                            Coexist26Bw(sp1,sp2)=sum(Bw(:,sp1).*Bw(:,sp2)>0)>0; %final coexistence matrix
                            Coexist26BLV1(sp1,sp2)=sum(BLV1(:,sp1).*BLV1(:,sp2)>0)>0; %final coexistence matrix
                            Coexist26BLV1w(sp1,sp2)=sum(BLV1w(:,sp1).*BLV1w(:,sp2)>0)>0; %final coexistence matrix
                        end
                    elseif (~isempty(find(Size5_6_Pos==sp1)))
                        if (~isempty(find(Size5_6_Pos==sp2))) %within size class coexistence
                            Coexist66Init(sp1,sp2)=sum(BtransEnd(:,sp1).*BtransEnd(:,sp2)>0)>0; %initial coexistence matrix
                            Coexist66B(sp1,sp2)=sum(B(:,sp1).*B(:,sp2)>0)>0; %final coexistence matrix
                            Coexist66Bw(sp1,sp2)=sum(Bw(:,sp1).*Bw(:,sp2)>0)>0; %final coexistence matrix
                            Coexist66BLV1(sp1,sp2)=sum(BLV1(:,sp1).*BLV1(:,sp2)>0)>0; %final coexistence matrix
                            Coexist66BLV1w(sp1,sp2)=sum(BLV1w(:,sp1).*BLV1w(:,sp2)>0)>0; %final coexistence matrix
                        end
                    end
                end
            end
            novelCoexistence=[novelCoexistence; sum((CoexistB(:)-CoexistInit(:))==1)/sum(CoexistB(:)) sum((CoexistBw(:)-CoexistInit(:))==1)/sum(CoexistBw(:)) sum((CoexistBLV1(:)-CoexistInit(:))==1)/sum(CoexistBLV1(:)) sum((CoexistBLV1w(:)-CoexistInit(:))==1)/sum(CoexistBLV1w(:))]; %portion of current coexisting pairs that were not there initially
            lostCoexistence=[lostCoexistence; sum((CoexistB(:)-CoexistInit(:))==-1)/sum(CoexistInit(:)) sum((CoexistBw(:)-CoexistInit(:))==-1)/sum(CoexistInit(:)) sum((CoexistBLV1(:)-CoexistInit(:))==-1)/sum(CoexistInit(:)) sum((CoexistBLV1w(:)-CoexistInit(:))==-1)/sum(CoexistInit(:))]; %portion of past coexisting pairs that are no longer there
            novelCoexistence22=[novelCoexistence22; sum((Coexist22B(:)-Coexist22Init(:))==1)/sum(Coexist22B(:)) sum((Coexist22Bw(:)-Coexist22Init(:))==1)/sum(Coexist22Bw(:)) sum((Coexist22BLV1(:)-Coexist22Init(:))==1)/sum(Coexist22BLV1(:)) sum((Coexist22BLV1w(:)-Coexist22Init(:))==1)/sum(Coexist22BLV1w(:))]; %portion of current coexisting pairs that were not there initially
            lostCoexistence22=[lostCoexistence22; sum((Coexist22B(:)-Coexist22Init(:))==-1)/sum(Coexist22Init(:)) sum((Coexist22Bw(:)-Coexist22Init(:))==-1)/sum(Coexist22Init(:)) sum((Coexist22BLV1(:)-Coexist22Init(:))==-1)/sum(Coexist22Init(:)) sum((Coexist22BLV1w(:)-Coexist22Init(:))==-1)/sum(Coexist22Init(:))]; %portion of past coexisting pairs that are no longer there
            novelCoexistence26=[novelCoexistence26; sum((Coexist26B(:)-Coexist26Init(:))==1)/sum(Coexist26B(:)) sum((Coexist26Bw(:)-Coexist26Init(:))==1)/sum(Coexist26Bw(:)) sum((Coexist26BLV1(:)-Coexist26Init(:))==1)/sum(Coexist26BLV1(:)) sum((Coexist26BLV1w(:)-Coexist26Init(:))==1)/sum(Coexist26BLV1w(:))]; %portion of current coexisting pairs that were not there initially
            lostCoexistence26=[lostCoexistence26; sum((Coexist26B(:)-Coexist26Init(:))==-1)/sum(Coexist26Init(:)) sum((Coexist26Bw(:)-Coexist26Init(:))==-1)/sum(Coexist26Init(:)) sum((Coexist26BLV1(:)-Coexist26Init(:))==-1)/sum(Coexist26Init(:)) sum((Coexist26BLV1w(:)-Coexist26Init(:))==-1)/sum(Coexist26Init(:))]; %portion of past coexisting pairs that are no longer there
            novelCoexistence66=[novelCoexistence66; sum((Coexist66B(:)-Coexist66Init(:))==1)/sum(Coexist66B(:)) sum((Coexist66Bw(:)-Coexist66Init(:))==1)/sum(Coexist66Bw(:)) sum((Coexist66BLV1(:)-Coexist66Init(:))==1)/sum(Coexist66BLV1(:)) sum((Coexist66BLV1w(:)-Coexist66Init(:))==1)/sum(Coexist66BLV1w(:))]; %portion of current coexisting pairs that were not there initially
            lostCoexistence66=[lostCoexistence66; sum((Coexist66B(:)-Coexist66Init(:))==-1)/sum(Coexist66Init(:)) sum((Coexist66Bw(:)-Coexist66Init(:))==-1)/sum(Coexist66Init(:)) sum((Coexist66BLV1(:)-Coexist66Init(:))==-1)/sum(Coexist66Init(:)) sum((Coexist66BLV1w(:)-Coexist66Init(:))==-1)/sum(Coexist66Init(:))]; %portion of past coexisting pairs that are no longer there
            
            
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
            meanInitCentroid=nansum(BtransEnd.*[1:numPatches]')./nansum(BtransEnd);
            meanInitCentroid(nansum(BtransEnd)==0)=NaN;
            InitCentroid=[InitCentroid; meanInitCentroid]; %mean location of each species at the end of transcient period
            meanFinalCentroid=nansum(Bnan.*[1:numPatches]')./nansum(Bnan);
            meanFinalCentroid(nansum(Bnan)==0)=NaN;
            meanFinalLVCentroid=nansum(BLV1nan.*[1:numPatches]')./nansum(BLV1nan);
            meanFinalLVCentroid(nansum(BLV1nan)==0)=NaN;
            FinalCentroid=[FinalCentroid; meanFinalCentroid]; %mean location of each species at the end of warming period
            FinalLVCentroid=[FinalLVCentroid; meanFinalLVCentroid]; %projected mean location of each species at the end of warming period
            meanFinalwCentroid=nansum(Bwnan.*[1:numPatches]')./nansum(Bwnan);
            meanFinalwCentroid(nansum(Bwnan)==0)=NaN;
            meanFinalwLVCentroid=nansum(BLV1wnan.*[1:numPatches]')./nansum(BLV1wnan);
            meanFinalwLVCentroid(nansum(BLV1wnan)==0)=NaN;
            FinalwCentroid=[FinalwCentroid; meanFinalwCentroid]; %mean location of each species at the end of warming period
            FinalwLVCentroid=[FinalwLVCentroid; meanFinalwLVCentroid]; %projected mean location of each species at the end of warming period
            CentroidShift0=[CentroidShift0; nanmean(meanFinalCentroid-meanInitCentroid) nanmean(meanFinalLVCentroid-meanInitCentroid)]; %average centroid shift
            CentroidShift=[CentroidShift; nanmean(meanFinalwCentroid-meanInitCentroid) nanmean(meanFinalwLVCentroid-meanInitCentroid)]; %average centroid shift
            Centroid2Shift0=[Centroid2Shift0; nanmean(meanFinalCentroid(Size2_3_Pos)-meanInitCentroid(Size2_3_Pos)) nanmean(meanFinalLVCentroid(Size2_3_Pos)-meanInitCentroid(Size2_3_Pos))]; %average centroid shift
            Centroid2Shift=[Centroid2Shift; nanmean(meanFinalwCentroid(Size2_3_Pos)-meanInitCentroid(Size2_3_Pos)) nanmean(meanFinalwLVCentroid(Size2_3_Pos)-meanInitCentroid(Size2_3_Pos))]; %average centroid shift
            Centroid6Shift0=[Centroid6Shift0; nanmean(meanFinalCentroid(Size5_6_Pos)-meanInitCentroid(Size5_6_Pos)) nanmean(meanFinalLVCentroid(Size5_6_Pos)-meanInitCentroid(Size5_6_Pos))]; %average centroid shift
            Centroid6Shift=[Centroid6Shift; nanmean(meanFinalwCentroid(Size5_6_Pos)-meanInitCentroid(Size5_6_Pos)) nanmean(meanFinalwLVCentroid(Size5_6_Pos)-meanInitCentroid(Size5_6_Pos))]; %average centroid shift
            All_CentroidShifts=cat(3,All_CentroidShifts,[meanFinalCentroid-meanInitCentroid;meanFinalwCentroid-meanInitCentroid;meanFinalLVCentroid-meanInitCentroid;meanFinalwLVCentroid-meanInitCentroid]);
            All_Centroids=cat(3,All_Centroids,meanInitCentroid);
            
            %trailing and leading edge shift
            leadingQ=0.025; %lower quantile indicates location of leading edge (towards cold region)
            trailingQ=0.975; %upper quantile indicates location of trailing edge (towards hot region)
            minSpeciesBiomass=leadingQ*nansum(BtransEnd); %determine leading biomass density
            maxSpeciesBiomass=trailingQ*nansum(BtransEnd); %determine trailing biomass density
            tempInitTrailing=sum(cumsum(BtransEnd,'omitnan')<maxSpeciesBiomass)+1;
            tempInitTrailing(nansum(BtransEnd)==0)=NaN;
            tempInitLeading=sum(cumsum(BtransEnd,'omitnan')<minSpeciesBiomass)+1;
            tempInitLeading(nansum(BtransEnd)==0)=NaN;
            InitTrailing=[InitTrailing; tempInitTrailing]; %mean location of each species at the end of transcient period
            InitLeading=[InitLeading; tempInitLeading]; %mean location of each species at the end of transcient period
            
            minSpeciesBiomass=leadingQ*nansum(Bnan); %determine trailing biomass density
            maxSpeciesBiomass=trailingQ*nansum(Bnan); %determine leading biomass density
            tempFinalTrailing=sum(cumsum(Bnan,'omitnan')<maxSpeciesBiomass)+1;
            tempFinalTrailing(nansum(Bnan)==0)=NaN;
            tempFinalLeading=sum(cumsum(Bnan,'omitnan')<minSpeciesBiomass)+1;
            tempFinalLeading(nansum(Bnan)==0)=NaN;
            FinalTrailing=[FinalTrailing; tempFinalTrailing]; %mean location of each species at the end of transcient period
            FinalLeading=[FinalLeading; tempFinalLeading]; %mean location of each species at the end of transcient period
            
            minSpeciesBiomass=leadingQ*nansum(Bwnan); %determine trailing biomass density
            maxSpeciesBiomass=trailingQ*nansum(Bwnan); %determine leading biomass density
            tempFinalwTrailing=sum(cumsum(Bwnan,'omitnan')<maxSpeciesBiomass)+1;
            tempFinalwTrailing(nansum(Bwnan)==0)=NaN;
            tempFinalwLeading=sum(cumsum(Bwnan,'omitnan')<minSpeciesBiomass)+1;
            tempFinalwLeading(nansum(Bwnan)==0)=NaN;
            FinalwTrailing=[FinalwTrailing; tempFinalwTrailing]; %mean location of each species at the end of transcient period
            FinalwLeading=[FinalwLeading; tempFinalwLeading]; %mean location of each species at the end of transcient period
            
            minSpeciesBiomass=leadingQ*nansum(BLV1nan); %determine trailing biomass density
            maxSpeciesBiomass=trailingQ*nansum(BLV1nan); %determine leading biomass density
            tempFinalLVTrailing=sum(cumsum(BLV1nan,'omitnan')<maxSpeciesBiomass)+1;
            tempFinalLVTrailing(nansum(BLV1nan)==0)=NaN;
            tempFinalLVLeading=sum(cumsum(BLV1nan,'omitnan')<minSpeciesBiomass)+1;
            tempFinalLVLeading(nansum(BLV1nan)==0)=NaN;
            FinalLVTrailing=[FinalLVTrailing; tempFinalLVTrailing]; %mean location of each species at the end of transcient period
            FinalLVLeading=[FinalLVLeading; tempFinalLVLeading]; %mean location of each species at the end of transcient period
            
            minSpeciesBiomass=leadingQ*nansum(BLV1wnan); %determine trailing biomass density
            maxSpeciesBiomass=trailingQ*nansum(BLV1wnan); %determine leading biomass density
            tempFinalwLVTrailing=sum(cumsum(BLV1wnan,'omitnan')<maxSpeciesBiomass)+1;
            tempFinalwLVTrailing(nansum(BLV1wnan)==0)=NaN;
            tempFinalwLVLeading=sum(cumsum(BLV1wnan,'omitnan')<minSpeciesBiomass)+1;
            tempFinalwLVLeading(nansum(BLV1wnan)==0)=NaN;
            FinalwLVTrailing=[FinalwLVTrailing; tempFinalwLVTrailing]; %mean location of each species at the end of transcient period
            FinalwLVLeading=[FinalwLVLeading; tempFinalwLVLeading]; %mean location of each species at the end of transcient period
            
            TrailingShift0=[TrailingShift0; nanmean(tempFinalTrailing-tempInitTrailing) nanmean(tempFinalLVTrailing-tempInitTrailing)]; %average trailing edge shift under no warming
            TrailingShift=[TrailingShift; nanmean(tempFinalwTrailing-tempInitTrailing) nanmean(tempFinalwLVTrailing-tempInitTrailing)]; %average trailing edge shift under warming
            LeadingShift0=[LeadingShift0; nanmean(tempFinalLeading-tempInitLeading) nanmean(tempFinalLVLeading-tempInitLeading)]; %average trailing edge shift under no warming
            LeadingShift=[LeadingShift; nanmean(tempFinalwLeading-tempInitLeading) nanmean(tempFinalwLVLeading-tempInitLeading)]; %average trailing edge shift under warming
            RangeExpansion0=[RangeExpansion0; nanmean((tempFinalTrailing-tempFinalLeading)-(tempInitTrailing-tempInitLeading)) nanmean((tempFinalLVTrailing-tempFinalLVLeading)-(tempInitTrailing-tempInitLeading))]; %average range expansion under no warming
            RangeExpansion=[RangeExpansion; nanmean((tempFinalwTrailing-tempFinalwLeading)-(tempInitTrailing-tempInitLeading)) nanmean((tempFinalwLVTrailing-tempFinalwLVLeading)-(tempInitTrailing-tempInitLeading))]; %average range expansion under warming
            Trailing2Shift0=[Trailing2Shift0; nanmean(tempFinalTrailing(Size2_3_Pos)-tempInitTrailing(Size2_3_Pos)) nanmean(tempFinalLVTrailing(Size2_3_Pos)-tempInitTrailing(Size2_3_Pos))]; %average trailing edge shift under no warming
            Trailing2Shift=[Trailing2Shift; nanmean(tempFinalwTrailing(Size2_3_Pos)-tempInitTrailing(Size2_3_Pos)) nanmean(tempFinalwLVTrailing(Size2_3_Pos)-tempInitTrailing(Size2_3_Pos))]; %average trailing edge shift under warming
            Leading2Shift0=[Leading2Shift0; nanmean(tempFinalLeading(Size2_3_Pos)-tempInitLeading(Size2_3_Pos)) nanmean(tempFinalLVLeading(Size2_3_Pos)-tempInitLeading(Size2_3_Pos))]; %average trailing edge shift under no warming
            Leading2Shift=[Leading2Shift; nanmean(tempFinalwLeading(Size2_3_Pos)-tempInitLeading(Size2_3_Pos)) nanmean(tempFinalwLVLeading(Size2_3_Pos)-tempInitLeading(Size2_3_Pos))]; %average trailing edge shift under warming
            Range2Expansion0=[Range2Expansion0; nanmean((tempFinalTrailing(Size2_3_Pos)-tempFinalLeading(Size2_3_Pos))-(tempInitTrailing(Size2_3_Pos)-tempInitLeading(Size2_3_Pos))) nanmean((tempFinalLVTrailing(Size2_3_Pos)-tempFinalLVLeading(Size2_3_Pos))-(tempInitTrailing(Size2_3_Pos)-tempInitLeading(Size2_3_Pos)))]; %average range expansion under no warming
            Range2Expansion=[Range2Expansion; nanmean((tempFinalwTrailing(Size2_3_Pos)-tempFinalwLeading(Size2_3_Pos))-(tempInitTrailing(Size2_3_Pos)-tempInitLeading(Size2_3_Pos))) nanmean((tempFinalwLVTrailing(Size2_3_Pos)-tempFinalwLVLeading(Size2_3_Pos))-(tempInitTrailing(Size2_3_Pos)-tempInitLeading(Size2_3_Pos)))]; %average range expansion under warming
            Trailing6Shift0=[Trailing6Shift0; nanmean(tempFinalTrailing(Size5_6_Pos)-tempInitTrailing(Size5_6_Pos)) nanmean(tempFinalLVTrailing(Size5_6_Pos)-tempInitTrailing(Size5_6_Pos))]; %average trailing edge shift under no warming
            Trailing6Shift=[Trailing6Shift; nanmean(tempFinalwTrailing(Size5_6_Pos)-tempInitTrailing(Size5_6_Pos)) nanmean(tempFinalwLVTrailing(Size5_6_Pos)-tempInitTrailing(Size5_6_Pos))]; %average trailing edge shift under warming
            Leading6Shift0=[Leading6Shift0; nanmean(tempFinalLeading(Size5_6_Pos)-tempInitLeading(Size5_6_Pos)) nanmean(tempFinalLVLeading(Size5_6_Pos)-tempInitLeading(Size5_6_Pos))]; %average trailing edge shift under no warming
            Leading6Shift=[Leading6Shift; nanmean(tempFinalwLeading(Size5_6_Pos)-tempInitLeading(Size5_6_Pos)) nanmean(tempFinalwLVLeading(Size5_6_Pos)-tempInitLeading(Size5_6_Pos))]; %average trailing edge shift under warming
            Range6Expansion0=[Range6Expansion0; nanmean((tempFinalTrailing(Size5_6_Pos)-tempFinalLeading(Size5_6_Pos))-(tempInitTrailing(Size5_6_Pos)-tempInitLeading(Size5_6_Pos))) nanmean((tempFinalLVTrailing(Size5_6_Pos)-tempFinalLVLeading(Size5_6_Pos))-(tempInitTrailing(Size5_6_Pos)-tempInitLeading(Size5_6_Pos)))]; %average range expansion under no warming
            Range6Expansion=[Range6Expansion; nanmean((tempFinalwTrailing(Size5_6_Pos)-tempFinalwLeading(Size5_6_Pos))-(tempInitTrailing(Size5_6_Pos)-tempInitLeading(Size5_6_Pos))) nanmean((tempFinalwLVTrailing(Size5_6_Pos)-tempFinalwLVLeading(Size5_6_Pos))-(tempInitTrailing(Size5_6_Pos)-tempInitLeading(Size5_6_Pos)))]; %average range expansion under warming
            
            RangeExpansionPerc0=[RangeExpansionPerc0; nanmean(((tempFinalTrailing-tempFinalLeading)-(tempInitTrailing-tempInitLeading))./(tempInitTrailing-tempInitLeading+1))*100 nanmean(((tempFinalLVTrailing-tempFinalLVLeading)-(tempInitTrailing-tempInitLeading))./(tempInitTrailing-tempInitLeading+1))*100]; %average range expansion % under no warming
            RangeExpansionPerc=[RangeExpansionPerc; nanmean(((tempFinalwTrailing-tempFinalwLeading)-(tempInitTrailing-tempInitLeading))./(tempInitTrailing-tempInitLeading+1))*100 nanmean(((tempFinalwLVTrailing-tempFinalwLVLeading)-(tempInitTrailing-tempInitLeading))./(tempInitTrailing-tempInitLeading+1))*100]; %average range expansion % under no warming
            Range2ExpansionPerc0=[Range2ExpansionPerc0; nanmean(((tempFinalTrailing(Size2_3_Pos)-tempFinalLeading(Size2_3_Pos))-(tempInitTrailing(Size2_3_Pos)-tempInitLeading(Size2_3_Pos)))./(tempInitTrailing(Size2_3_Pos)-tempInitLeading(Size2_3_Pos)+1))*100 nanmean(((tempFinalLVTrailing(Size2_3_Pos)-tempFinalLVLeading(Size2_3_Pos))-(tempInitTrailing(Size2_3_Pos)-tempInitLeading(Size2_3_Pos)))./(tempInitTrailing(Size2_3_Pos)-tempInitLeading(Size2_3_Pos)+1))*100]; %average range expansion % under no warming
            Range2ExpansionPerc=[Range2ExpansionPerc; nanmean(((tempFinalwTrailing(Size2_3_Pos)-tempFinalwLeading(Size2_3_Pos))-(tempInitTrailing(Size2_3_Pos)-tempInitLeading(Size2_3_Pos)))./(tempInitTrailing(Size2_3_Pos)-tempInitLeading(Size2_3_Pos)+1))*100 nanmean(((tempFinalwLVTrailing(Size2_3_Pos)-tempFinalwLVLeading(Size2_3_Pos))-(tempInitTrailing(Size2_3_Pos)-tempInitLeading(Size2_3_Pos)))./(tempInitTrailing(Size2_3_Pos)-tempInitLeading(Size2_3_Pos)+1))*100]; %average range expansion % under no warming
            Range6ExpansionPerc0=[Range6ExpansionPerc0; nanmean(((tempFinalTrailing(Size5_6_Pos)-tempFinalLeading(Size5_6_Pos))-(tempInitTrailing(Size5_6_Pos)-tempInitLeading(Size5_6_Pos)))./(tempInitTrailing(Size5_6_Pos)-tempInitLeading(Size5_6_Pos)+1))*100 nanmean(((tempFinalLVTrailing(Size5_6_Pos)-tempFinalLVLeading(Size5_6_Pos))-(tempInitTrailing(Size5_6_Pos)-tempInitLeading(Size5_6_Pos)))./(tempInitTrailing(Size5_6_Pos)-tempInitLeading(Size5_6_Pos)+1))*100]; %average range expansion % under no warming
            Range6ExpansionPerc=[Range6ExpansionPerc; nanmean(((tempFinalwTrailing(Size5_6_Pos)-tempFinalwLeading(Size5_6_Pos))-(tempInitTrailing(Size5_6_Pos)-tempInitLeading(Size5_6_Pos)))./(tempInitTrailing(Size5_6_Pos)-tempInitLeading(Size5_6_Pos)+1))*100 nanmean(((tempFinalwLVTrailing(Size5_6_Pos)-tempFinalwLVLeading(Size5_6_Pos))-(tempInitTrailing(Size5_6_Pos)-tempInitLeading(Size5_6_Pos)))./(tempInitTrailing(Size5_6_Pos)-tempInitLeading(Size5_6_Pos)+1))*100]; %average range expansion % under no warming
            
            All_TrailingShifts=cat(3,All_TrailingShifts,[tempFinalTrailing-tempInitTrailing;tempFinalwTrailing-tempInitTrailing;tempFinalLVTrailing-tempInitTrailing;tempFinalwLVTrailing-tempInitTrailing]);
            All_LeadingShifts=cat(3,All_LeadingShifts,[tempFinalLeading-tempInitLeading;tempFinalwLeading-tempInitLeading;tempFinalLVLeading-tempInitLeading;tempFinalwLVLeading-tempInitLeading]);
            All_RangeExpansions=cat(3,All_RangeExpansions,[(tempFinalTrailing-tempFinalLeading)-(tempInitTrailing-tempInitLeading);(tempFinalwTrailing-tempFinalwLeading)-(tempInitTrailing-tempInitLeading);(tempFinalLVTrailing-tempFinalLVLeading)-(tempInitTrailing-tempInitLeading);(tempFinalwLVTrailing-tempFinalwLVLeading)-(tempInitTrailing-tempInitLeading)]);
            All_RangeSizes=cat(3,All_RangeSizes,tempInitTrailing-tempInitLeading+1);
        end
    end
end

%FileName=sprintf('WarmingMovementStats_%s.mat', TimeData);
FileName=sprintf('WarmingDispersalStats.mat');
save(FileName);
