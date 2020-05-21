function [LL,BiomR2,ProdR2,r1,a1,z1,BLV1,gainBLV1,BLV1w,gainBLV1w]=Make_warming_SingleSpeciesProjection_LL(fitCode,Btrans,dBtrans,gainBtrans,B_yrs,gainB_yrs,P)

%numIt=40;    %number of iterations for each parameter combination
SampInt=365;  %record every SampInt pts (days) in time series during transcient period

TimeData=string(datetime);

%timepoints to record:
RecordYrStartRange=[1600 2400] %warming starts randomly between the two years 1600-2400
TransEndYrs=800; %number of years at the end of transcient (burn-in) period used to fit 800
make_traits;
numT=length(tempChange); %number of temperature change scenarios

%RecordYrStart=round((RecordYrStartRange(2)-RecordYrStartRange(1))*rand+RecordYrStartRange(1)); %randomize initial time of warming (rounded to day)
RecordYrStart=(P.Tend-1)/365-200; %infer from recorded end year
TimePts=[RecordYrStart*365+1:365:(RecordYrStart+200)*365+1]; %record every year
numPts=length(TimePts); %number of time points

%    [P B Z T] = make_parameters(TR,i);
P.Tend=TimePts(end); %set last point to be 200 years after end of random-length transcient period (in days)
temps=num2str(P.dT*73000);
sprintf('Foodweb_%s_%s%s.mat', TimeData, num2str(P.iter), ['_numSpecies' num2str(P.n) '_dT' temps(~isspace(temps)) '_basalSize' num2str(P.s.m0) '_meanD' num2str(HP.sdm(i)) '_stdD' num2str(HP.sdv(i))])


r1=0;
a1=0;
z1=0;
T1w=[];
BLV1=[];
BLV1w=[];
gainBLV1w=zeros(P.nx,P.n,numT);
dBLV1w=zeros(P.nx,P.n,numT);

BLV1_yrs=zeros(P.nx,P.n,numPts);
gainBLV1_yrs=zeros(P.nx,P.n,numPts);
BLV1w_yrs=zeros(P.nx,P.n,numPts,numT);
gainBLV1w_yrs=zeros(P.nx,P.n,numPts,numT);

%fitCode=[.5 .5];
[r1,a1,z1,K1,flag1,raR21,r_T1,K_T1,K_T_ratio1,r_T_ratio1]=estSingleSpeciesModelmsy(Btrans,dBtrans,gainBtrans,P,fitCode(1,:)); %fit growth model to all patches at once
%
%% Iterate model forward in time (days at the moment)
% YearStartT=1;
% ttrans=1;
% for t = 1:P.Tend
%     if t==TimePts(1) %transition time when no-warming and warming experiments diverge from common states
%         BLV1=repmat(Btrans(:,:,end),1,1); %no warming case under estimated single-species dynamics
%         BLV1w=repmat(Btrans(:,:,end),1,1,numT); %warming case under estimated single-species dynamics
%     end
%     
%     if t>=TimePts(1)
%         BLV1(BLV1<eps) = 0;
%         BLV1w(BLV1w<eps)= 0;
%     end
%     
%     
%     % Shift thermal gradient
%     T1      = P.T;% + t.*P.dT; %<<< add this when time is right
%     if t>=TimePts(1)
%         T1w      = repmat(P.T,1,numT) + (t-TimePts(1)-1)*P.dT;
%     end
%     %T1      = P.T + 2*sin(2*pi*((365-t.*P.dt)/365)); %seasonal cycle
%     %T1      = P.T + 2*sin(2*pi*((365-t.*P.dt/10)/365)); %decadal cycle
%     
%     if t>=TimePts(1)
%         BLV1=sub_move(BLV1,P); %single species model no temp change
%         [gainBLV1 dBLV1] = sub_demogLVmsy(BLV1,T1,r1,a1,z1,P);
%         for TCase=1:numT %for each temperature change scenario
%             BLV1w(:,:,TCase)        = sub_move(BLV1w(:,:,TCase),P); % move
%             [gainBLV1w(:,:,TCase) dBLV1w(:,:,TCase)] = sub_demogLVmsy(BLV1w(:,:,TCase),T1w(:,TCase),r1,a1,z1,P); % grow/die
%         end
%     end
%     
%     %record major time points:
%     tpos=find(t==TimePts);
%     if ~isempty(tpos)
%         BLV1_yrs(:,:,tpos)=BLV1;
%         gainBLV1_yrs(:,:,tpos)=gainBLV1;
%         BLV1w_yrs(:,:,tpos,:)=BLV1w;
%         gainBLV1w_yrs(:,:,tpos,:)=gainBLV1w;
%     end
%     
%     if t>=TimePts(1)
%         BLV1 = BLV1 + dBLV1;
%         BLV1w = BLV1w + dBLV1w;
%     end
% end

%     FoodWebFile=sprintf('Foodweb_%s_%s%s.mat', TimeData, num2str(P.iter), ['_numSpecies' num2str(P.n) '_dT' num2str(P.dT*73000) '_basalSize' num2str(P.s.m0) '_meanD' num2str(HP.sdm(i)) '_stdD' num2str(HP.sdv(i))]);
%     %savemat_Foodweb(FoodWebFile,Z_yrs,B_yrs,gainB_yrs,gainZ_yrs,v_yrs,TE_yrs,PB_yrs,TLik_yrs,TLi_yrs,TLk_yrs,TLall_yrs,Zw_yrs,Bw_yrs,gainBw_yrs,gainZw_yrs,vw_yrs,TEw_yrs,PBw_yrs,TLikw_yrs,TLiw_yrs,TLkw_yrs,TLallw_yrs,BLV_yrs,gainBLV_yrs,BLVw_yrs,gainBLVw_yrs,BLV1_yrs,gainBLV1_yrs,BLV1w_yrs,gainBLV1w_yrs,BLV2_yrs,gainBLV2_yrs,BLV2w_yrs,gainBLV2w_yrs,BLV3_yrs,gainBLV3_yrs,BLV3w_yrs,gainBLV3w_yrs,BLV4_yrs,gainBLV4_yrs,BLV4w_yrs,gainBLV4w_yrs,BLV5_yrs,gainBLV5_yrs,BLV5w_yrs,gainBLV5w_yrs,BLV6_yrs,gainBLV6_yrs,BLV6w_yrs,gainBLV6w_yrs,Btrans,dBtrans,gainBtrans,r,a,z,K,flag,raR2,r_T,K_T,K_T_ratio,r_T_ratio,r1,a1,z1,K1,flag1,raR21,r_T1,K_T1,K_T_ratio1,r_T_ratio1,r2,a2,z2,K2,flag2,raR22,r_T2,K_T2,K_T_ratio2,r_T_ratio2,r3,a3,z3,K3,flag3,raR23,r_T3,K_T3,K_T_ratio3,r_T_ratio3,r4,a4,z4,K4,flag4,raR24,r_T4,K_T4,K_T_ratio4,r_T_ratio4,r5,a5,z5,K5,flag5,raR25,r_T5,K_T5,K_T_ratio5,r_T_ratio5,r6,a6,z6,K6,flag6,raR26,r_T6,K_T6,K_T_ratio6,r_T_ratio6,fitCode,P);
%     savemat_Foodweb(FoodWebFile,Z_yrs,B_yrs,gainB_yrs,gainZ_yrs,v_yrs,TE_yrs,PB_yrs,TLik_yrs,TLi_yrs,TLk_yrs,TLall_yrs,Zw_yrs,Bw_yrs,gainBw_yrs,gainZw_yrs,vw_yrs,TEw_yrs,PBw_yrs,TLikw_yrs,TLiw_yrs,TLkw_yrs,TLallw_yrs,BLV1_yrs,gainBLV1_yrs,BLV1w_yrs,gainBLV1w_yrs,BLV2_yrs,gainBLV2_yrs,BLV2w_yrs,gainBLV2w_yrs,BLV3_yrs,gainBLV3_yrs,BLV3w_yrs,gainBLV3w_yrs,BLV4_yrs,gainBLV4_yrs,BLV4w_yrs,gainBLV4w_yrs,Btrans,dBtrans,gainBtrans,r1,a1,z1,K1,flag1,raR21,r_T1,K_T1,K_T_ratio1,r_T_ratio1,r2,a2,z2,K2,flag2,raR22,r_T2,K_T2,K_T_ratio2,r_T_ratio2,r3,a3,z3,K3,flag3,raR23,r_T3,K_T3,K_T_ratio3,r_T_ratio3,r4,a4,z4,K4,flag4,raR24,r_T4,K_T4,K_T_ratio4,r_T_ratio4,fitCode,P);
