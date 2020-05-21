%Make_warming_MultPtsstats_Parallel.m
%Edward Tekwa Oct 22, 17
%run food web simulations with parallel warming and no-warming cases
%Nov 3, 17: run multiple temperature change scenarios together and record
%several time points for each simulation, on parallel clusters

% %clear all; %close all;
delete(gcp('nocreate'))
% mycluster=parcluster;
% parpool(min(mycluster.NumWorkers,10)) %run on either the max number of clusters or the limit specified here
parpool('local', str2num(getenv('SLURM_CPUS_PER_TASK')))
warning('off','all')

numIt=40;    %number of iterations for each parameter combination
SampInt=365;  %record every SampInt pts (days) in time series during transcient period

TimeData=string(datetime);

%timepoints to record:
%TimePts=[100000 118250 136500 154750 173000]; %first point is when temperature change starts; last point is end of simulation (still warming)
%TimePts=[219000:365:292000]; %record every year (600-800 yrs)
%RecordYrStart=2000;
RecordYrStartRange=[1600 2400] %warming starts randomly between the two years 1600-2400
% RecordYrStart=
% TimePts=[RecordYrStart*365+1:365:(RecordYrStart+200)*365+1]; %record every year
%TimePts=[365000:365:438000]; %record every year (1000-1200 yrs)
TransEndYrs=800; %number of years at the end of transcient (burn-in) period used to fit 800
%single-species model
% numPts=length(TimePts); %number of time points
make_traits;
numT=length(tempChange); %number of temperature change scenarios

numSims=numIt*length(specialist)*length(sdm)
estSimDays=(numSims/str2num(getenv('SLURM_CPUS_PER_TASK')))*12/24 %takes about 12 hrs for 1 simulation per core on Amarel

%parfor
parfor i = 1:size(TR,2) %go through all parameter combinations (as defined in make_traits.m)
    %% make parameters for a given set of traits
    RecordYrStart=round((RecordYrStartRange(2)-RecordYrStartRange(1))*rand+RecordYrStartRange(1)); %randomize initial time of warming (rounded to day)
    TimePts=[RecordYrStart*365+1:365:(RecordYrStart+200)*365+1]; %record every year
    numPts=length(TimePts); %number of time points

    [P B Z T] = make_parameters(TR,i);
    P.Tend=TimePts(end); %set last point to be 200 years after end of random-length transcient period (in days)
    temps=num2str(P.dT*73000);
    sprintf('Foodweb_%s_%s%s.mat', TimeData, num2str(P.iter), ['_numSpecies' num2str(P.n) '_dT' temps(~isspace(temps)) '_basalSize' num2str(P.s.m0) '_meanD' num2str(HP.sdm(i)) '_stdD' num2str(HP.sdv(i))])
    
%     r=0;
%     rs=0;
%     a=0;
%     z=0;
    r1=0;
    a1=0;
    z1=0;
    r2=0;
    a2=0;
    z2=0;
    r3=0;
    a3=0;
    z3=0;
    r4=0;
    a4=0;
    z4=0;
%     r5=0;
%     a5=0;
%     z5=0;
%     r6=0;
%     a6=0;
%     z6=0;
    K=0;
    flag=0;
    raR2=0;
    r_T=0;
    K_T_ratio=0;
    r_T_ratio=0;
    
    T1w=[];
    Bw=[];
    Zw=[];
%     BLV=[];
%     BLVw=[];
    BLV1=[];
    BLV1w=[];
    BLV2=[];
    BLV2w=[];
    BLV3=[];
    BLV3w=[];
    BLV4=[];
    BLV4w=[];
%     BLV5=[];
%     BLV5w=[];
%     BLV6=[];
%     BLV6w=[];
    gainBw=zeros(P.nx,P.n,numT);
    gainZw=zeros(P.nx,1,numT);
    dBw=zeros(P.nx,P.n,numT);
    dZw=zeros(P.nx,1,numT);
    vw=zeros(P.n, P.nx,numT);
    TEw=zeros(P.n, P.nx,numT);
    PBw=zeros(P.n, P.nx,numT);
    TLikw=zeros(P.nx,P.n+1,numT);
    TLiw=zeros(1,P.n,numT);
    TLkw=zeros(P.nx,1,numT);
    TLallw=zeros(1,numT);
%     gainBLVw=zeros(P.nx,P.n,numT);
%     dBLVw=zeros(P.nx,P.n,numT);
    gainBLV1w=zeros(P.nx,P.n,numT);
    dBLV1w=zeros(P.nx,P.n,numT);
    gainBLV2w=zeros(P.nx,P.n,numT);
    dBLV2w=zeros(P.nx,P.n,numT);
    gainBLV3w=zeros(P.nx,P.n,numT);
    dBLV3w=zeros(P.nx,P.n,numT);
    gainBLV4w=zeros(P.nx,P.n,numT);
    dBLV4w=zeros(P.nx,P.n,numT);
%     gainBLV5w=zeros(P.nx,P.n,numT);
%     dBLV5w=zeros(P.nx,P.n,numT);
%     gainBLV6w=zeros(P.nx,P.n,numT);
%     dBLV6w=zeros(P.nx,P.n,numT);
    
    Z_yrs=zeros(P.nx,1,numPts);
    B_yrs=zeros(P.nx,P.n,numPts);
    TE_yrs = zeros(P.n, P.nx,numPts); %trophic efficiency
    PB_yrs = zeros(P.n, P.nx,numPts); %doubling time
    gainZ_yrs=zeros(P.nx, 1,numPts); %productivity of basal resource
    gainB_yrs=zeros(P.nx, P.n,numPts); %productivity of heterotrophs
    v_yrs=zeros(P.n,P.nx,numPts);
    TLik_yrs=zeros(P.nx,P.n+1,numPts);
    TLi_yrs=zeros(1,P.n,numPts);
    TLk_yrs=zeros(P.nx,1,numPts);
    TLall_yrs=zeros(numPts,1);
%     BLV_yrs=zeros(P.nx,P.n,numPts);
%     gainBLV_yrs=zeros(P.nx,P.n,numPts);
    BLV1_yrs=zeros(P.nx,P.n,numPts);
    gainBLV1_yrs=zeros(P.nx,P.n,numPts);
    BLV2_yrs=zeros(P.nx,P.n,numPts);
    gainBLV2_yrs=zeros(P.nx,P.n,numPts);
    BLV3_yrs=zeros(P.nx,P.n,numPts);
    gainBLV3_yrs=zeros(P.nx,P.n,numPts);
    BLV4_yrs=zeros(P.nx,P.n,numPts);
    gainBLV4_yrs=zeros(P.nx,P.n,numPts);
%     BLV5_yrs=zeros(P.nx,P.n,numPts);
%     gainBLV5_yrs=zeros(P.nx,P.n,numPts);
%     BLV6_yrs=zeros(P.nx,P.n,numPts);
%     gainBLV6_yrs=zeros(P.nx,P.n,numPts);
    
    %with warming: various degrees increase
    Zw_yrs=zeros(P.nx,1,numPts,numT);
    Bw_yrs=zeros(P.nx,P.n,numPts,numT);
    TEw_yrs = zeros(P.n, P.nx,numPts,numT);
    PBw_yrs = zeros(P.n, P.nx,numPts,numT);
    gainZw_yrs=zeros(P.nx, 1,numPts,numT);
    gainBw_yrs=zeros(P.nx, P.n,numPts,numT);
    vw_yrs=zeros(P.n,P.nx,numPts,numT);
    TLikw_yrs=zeros(P.nx,P.n+1,numPts,numT);
    TLiw_yrs=zeros(1,P.n,numPts,numT);
    TLkw_yrs=zeros(P.nx,1,numPts,numT);
    TLallw_yrs=zeros(numPts,numT);
%     BLVw_yrs=zeros(P.nx,P.n,numPts,numT);
%     gainBLVw_yrs=zeros(P.nx,P.n,numPts,numT);
    BLV1w_yrs=zeros(P.nx,P.n,numPts,numT);
    gainBLV1w_yrs=zeros(P.nx,P.n,numPts,numT);
    BLV2w_yrs=zeros(P.nx,P.n,numPts,numT);
    gainBLV2w_yrs=zeros(P.nx,P.n,numPts,numT);
    BLV3w_yrs=zeros(P.nx,P.n,numPts,numT);
    gainBLV3w_yrs=zeros(P.nx,P.n,numPts,numT);
    BLV4w_yrs=zeros(P.nx,P.n,numPts,numT);
    gainBLV4w_yrs=zeros(P.nx,P.n,numPts,numT);
%     BLV5w_yrs=zeros(P.nx,P.n,numPts,numT);
%     gainBLV5w_yrs=zeros(P.nx,P.n,numPts,numT);
%     BLV6w_yrs=zeros(P.nx,P.n,numPts,numT);
%     gainBLV6w_yrs=zeros(P.nx,P.n,numPts,numT);
    
    %matrix to record biomasses and changes in species biomasses per patch during transcient period
%     dBtrans=zeros(P.nx,P.n,(TimePts(1)-1)/(2*SampInt));
%     Btrans=zeros(P.nx,P.n,(TimePts(1)-1)/(2*SampInt));
%     gainBtrans=zeros(P.nx,P.n,(TimePts(1)-1)/(2*SampInt));
    dBtrans=zeros(P.nx,P.n,TransEndYrs);
    Btrans=zeros(P.nx,P.n,TransEndYrs);
    gainBtrans=zeros(P.nx,P.n,TransEndYrs);
    %% Iterate model forward in time (days at the moment)
    YearStartT=1;
    ttrans=1;
    for t = 1:P.Tend
        if t==TimePts(1) %transition time when no-warming and warming experiments diverge from common states
            %Bw=B; Zw =Z;
            Bw=repmat(B,1,1,numT);
            Zw=repmat(Z,1,1,numT);
%             BLV=repmat(B,1,1); %no warming case under estimated Lotka-Volterra dynamics
%             BLVw=repmat(B,1,1,numT); %warming case under estimated Lotka-Volterra dynamics
            BLV1=repmat(B,1,1); %no warming case under estimated single-species dynamics
            BLV1w=repmat(B,1,1,numT); %warming case under estimated single-species dynamics
            BLV2=repmat(B,1,1); %no warming case under estimated single-species dynamics
            BLV2w=repmat(B,1,1,numT); %warming case under estimated single-species dynamics
            BLV3=repmat(B,1,1); %no warming case under estimated single-species dynamics
            BLV3w=repmat(B,1,1,numT); %warming case under estimated single-species dynamics
            BLV4=repmat(B,1,1); %no warming case under estimated single-species dynamics
            BLV4w=repmat(B,1,1,numT); %warming case under estimated single-species dynamics
%             BLV5=repmat(B,1,1); %no warming case under estimated single-species dynamics
%             BLV5w=repmat(B,1,1,numT); %warming case under estimated single-species dynamics
%             BLV6=repmat(B,1,1); %no warming case under estimated single-species dynamics
%             BLV6w=repmat(B,1,1,numT); %warming case under estimated single-species dynamics
        end
        
        % fix
        B(B<eps) = 0;% eps or 0;
        Z(Z<eps) = 0;% eps;
        if t>=TimePts(1)
            Bw(Bw<eps) = 0;% eps;
            Zw(Zw<eps) = 0;% eps;
%             BLV(BLV<eps) = 0;
%             BLVw(BLVw<eps)= 0;
            BLV1(BLV1<eps) = 0;
            BLV1w(BLV1w<eps)= 0;
            BLV2(BLV2<eps) = 0;
            BLV2w(BLV2w<eps)= 0;
            BLV3(BLV3<eps) = 0;
            BLV3w(BLV3w<eps)= 0;
            BLV4(BLV1<eps) = 0;
            BLV4w(BLV1w<eps)= 0;
%             BLV5(BLV2<eps) = 0;
%             BLV5w(BLV2w<eps)= 0;
%             BLV6(BLV3<eps) = 0;
%             BLV6w(BLV3w<eps)= 0;
        end
        
        
        % Shift thermal gradient
        T1      = P.T;% + t.*P.dT; %<<< add this when time is right
        if t>=TimePts(1)
            T1w      = repmat(P.T,1,numT) + (t-TimePts(1)-1)*P.dT;
        end
        %T1      = P.T + 2*sin(2*pi*((365-t.*P.dt)/365)); %seasonal cycle
        %T1      = P.T + 2*sin(2*pi*((365-t.*P.dt/10)/365)); %decadal cycle
        
        
        % demographics
        B        = sub_move(B,P); % move
        [gainB gainZ dB dZ v TE PB TLik TLi TLk TLall] = sub_demog(t,B,Z,T1,P); % grow/die
        if t>=TimePts(1)
%            BLV=sub_move(BLV,P); %LV no temp change
            BLV1=sub_move(BLV1,P); %single species model no temp change
            BLV2=sub_move(BLV2,P); %single species model no temp change
            BLV3=sub_move(BLV3,P); %single species model no temp change
            BLV4=sub_move(BLV4,P); %single species model no temp change
%             BLV5=sub_move(BLV5,P); %single species model no temp change
%             BLV6=sub_move(BLV6,P); %single species model no temp change
%             [gainBLV dBLV] = sub_demogLV(BLV,T1,r,a,z,P); %grow/die according to Lotka-Volterra approximation
            [gainBLV1 dBLV1] = sub_demogLVmsy(BLV1,T1,r1,a1,z1,P);
            [gainBLV2 dBLV2] = sub_demogLVmsy(BLV2,T1,r2,a2,z2,P);
            [gainBLV3 dBLV3] = sub_demogLVmsy(BLV3,T1,r3,a3,z3,P);
            [gainBLV4 dBLV4] = sub_demogLVmsy(BLV4,T1,r4,a4,z4,P);
%             [gainBLV5 dBLV5] = sub_demogLVmsy(BLV5,T1,r5,a5,z5,P);
%             [gainBLV6 dBLV6] = sub_demogLVmsy(BLV6,T1,r6,a6,z6,P);
            for TCase=1:numT %for each temperature change scenario
                Bw(:,:,TCase)        = sub_move(Bw(:,:,TCase),P); % move
                [gainBw(:,:,TCase) gainZw(:,:,TCase) dBw(:,:,TCase) dZw(:,:,TCase) vw(:,:,TCase) TEw(:,:,TCase) PBw(:,:,TCase) TLikw(:,:,TCase) TLiw(:,:,TCase) TLkw(:,:,TCase) TLallw(TCase)] = sub_demog(t,Bw(:,:,TCase),Zw(:,:,TCase),T1w(:,TCase),P); % grow/die
%                 BLVw(:,:,TCase)        = sub_move(BLVw(:,:,TCase),P); % move
%                 [gainBLVw(:,:,TCase) dBLVw(:,:,TCase)] = sub_demogLV(BLVw(:,:,TCase),T1w(:,TCase),r,a,z,P); % grow/die
                BLV1w(:,:,TCase)        = sub_move(BLV1w(:,:,TCase),P); % move
                [gainBLV1w(:,:,TCase) dBLV1w(:,:,TCase)] = sub_demogLVmsy(BLV1w(:,:,TCase),T1w(:,TCase),r1,a1,z1,P); % grow/die
                BLV2w(:,:,TCase)        = sub_move(BLV2w(:,:,TCase),P); % move
                [gainBLV2w(:,:,TCase) dBLV2w(:,:,TCase)] = sub_demogLVmsy(BLV2w(:,:,TCase),T1w(:,TCase),r2,a2,z2,P); % grow/die
                BLV3w(:,:,TCase)        = sub_move(BLV3w(:,:,TCase),P); % move
                [gainBLV3w(:,:,TCase) dBLV3w(:,:,TCase)] = sub_demogLVmsy(BLV3w(:,:,TCase),T1w(:,TCase),r3,a3,z3,P); % grow/die
                BLV4w(:,:,TCase)        = sub_move(BLV4w(:,:,TCase),P); % move
                [gainBLV4w(:,:,TCase) dBLV4w(:,:,TCase)] = sub_demogLVmsy(BLV4w(:,:,TCase),T1w(:,TCase),r4,a4,z4,P); % grow/die
%                 BLV5w(:,:,TCase)        = sub_move(BLV5w(:,:,TCase),P); % move
%                 [gainBLV5w(:,:,TCase) dBLV5w(:,:,TCase)] = sub_demogLVmsy(BLV5w(:,:,TCase),T1w(:,TCase),r5,a5,z5,P); % grow/die
%                 BLV6w(:,:,TCase)        = sub_move(BLV6w(:,:,TCase),P); % move
%                 [gainBLV6w(:,:,TCase) dBLV6w(:,:,TCase)] = sub_demogLVmsy(BLV6w(:,:,TCase),T1w(:,TCase),r6,a6,z6,P); % grow/die

            end
        end
        
        %record transcient time points (every SampInt days):
        if t>(TimePts(1)-TransEndYrs*365) && t<=TimePts(1) &&  floor(t/SampInt)==t/SampInt %last condition is satisfied every SampInt points
        %if t>TimePts(1)/2 && t<=TimePts(1) &&  floor(t/SampInt)==t/SampInt %last condition is satisfied every SampInt points
            Btrans(:,:,ttrans)=B; %biomass (before current time step update)
            dBtrans(:,:,ttrans)=dB * P.dt; %net change in biomass
            gainBtrans(:,:,ttrans)=gainB; %productivity
            ttrans=ttrans+1;
        end
        
        %record major time points:
        tpos=find(t==TimePts);
        if ~isempty(tpos)
            Z_yrs(:,:,tpos)=Z;
            B_yrs(:,:,tpos)=B;
            Zw_yrs(:,:,tpos,:)=Zw;
            Bw_yrs(:,:,tpos,:)=Bw;
            gainB_yrs(:,:,tpos)=gainB;
            gainZ_yrs(:,:,tpos)=gainZ;
            v_yrs(:,:,tpos)=v;
            TE_yrs(:,:,tpos)=TE;
            PB_yrs(:,:,tpos)=PB;
            TLik_yrs(:,:,tpos)=TLik;
            TLi_yrs(:,:,tpos)=TLi;
            TLk_yrs(:,:,tpos)=TLk;
            TLall_yrs(tpos)=TLall;
            gainBw_yrs(:,:,tpos,:)=gainBw;
            gainZw_yrs(:,:,tpos,:)=gainZw;
            vw_yrs(:,:,tpos,:)=vw;
            TEw_yrs(:,:,tpos,:)=TEw;
            PBw_yrs(:,:,tpos,:)=PBw;
            TLikw_yrs(:,:,tpos,:)=TLikw;
            TLiw_yrs(:,:,tpos,:)=TLiw;
            TLkw_yrs(:,:,tpos,:)=TLkw;
            TLallw_yrs(tpos,:)=TLallw;
%             BLV_yrs(:,:,tpos)=BLV;
%             gainBLV_yrs(:,:,tpos)=gainBLV;
%             BLVw_yrs(:,:,tpos,:)=BLVw;
%             gainBLVw_yrs(:,:,tpos,:)=gainBLVw;
            BLV1_yrs(:,:,tpos)=BLV1;
            gainBLV1_yrs(:,:,tpos)=gainBLV1;
            BLV1w_yrs(:,:,tpos,:)=BLV1w;
            gainBLV1w_yrs(:,:,tpos,:)=gainBLV1w;
            BLV2_yrs(:,:,tpos)=BLV2;
            gainBLV2_yrs(:,:,tpos)=gainBLV2;
            BLV2w_yrs(:,:,tpos,:)=BLV2w;
            gainBLV2w_yrs(:,:,tpos,:)=gainBLV2w;
            BLV3_yrs(:,:,tpos)=BLV3;
            gainBLV3_yrs(:,:,tpos)=gainBLV3;
            BLV3w_yrs(:,:,tpos,:)=BLV3w;
            gainBLV3w_yrs(:,:,tpos,:)=gainBLV3w;
            BLV4_yrs(:,:,tpos)=BLV4;
            gainBLV4_yrs(:,:,tpos)=gainBLV4;
            BLV4w_yrs(:,:,tpos,:)=BLV4w;
            gainBLV4w_yrs(:,:,tpos,:)=gainBLV4w;
%             BLV5_yrs(:,:,tpos)=BLV5;
%             gainBLV5_yrs(:,:,tpos)=gainBLV5;
%             BLV5w_yrs(:,:,tpos,:)=BLV5w;
%             gainBLV5w_yrs(:,:,tpos,:)=gainBLV5w;
%             BLV6_yrs(:,:,tpos)=BLV6;
%             gainBLV6_yrs(:,:,tpos)=gainBLV6;
%             BLV6w_yrs(:,:,tpos,:)=BLV6w;
%             gainBLV6w_yrs(:,:,tpos,:)=gainBLV6w;
        end
        
        % update biomasses
        Z = Z + (dZ * P.dt);
        B = B + (dB * P.dt);
        if t>=TimePts(1)
            Zw = Zw + (dZw * P.dt);
            Bw = Bw + (dBw * P.dt);
%             BLV = BLV + dBLV;
%             BLVw = BLVw + dBLVw;
            BLV1 = BLV1 + dBLV1;
            BLV1w = BLV1w + dBLV1w;
            BLV2 = BLV2 + dBLV2;
            BLV2w = BLV2w + dBLV2w;
            BLV3 = BLV3 + dBLV3;
            BLV3w = BLV3w + dBLV3w;
            BLV4 = BLV4 + dBLV4;
            BLV4w = BLV4w + dBLV4w;
%             BLV5 = BLV5 + dBLV5;
%             BLV5w = BLV5w + dBLV5w;
%             BLV6 = BLV6 + dBLV6;
%             BLV6w = BLV6w + dBLV6w;
        end
        
        %estimate ri(T) and ki(T) for each species at different
        %temperatures, to be used in projection with no species
        %interactions
        if t==TimePts(1)-1
            %estimate per-mass intrinsic growth rate and competition terms
            %[r,a,K,flag,raR2,r_T,rs,K_T_ratio,r_T_ratio]=estSpeciesModel(Btrans,dBtrans,gainBtrans,P); %fit growth model to each patch independently, then estimate intrinsic growth as function of T
            %[r1,a1,K1,flag1,raR21,r_T1,K_T1,K_T_ratio1,r_T_ratio1]=estSingleSpeciesModel(Btrans(end-TransEndYrs:end),dBtrans(end-TransEndYrs:end),gainBtrans(end-TransEndYrs:end),P); %fit growth model to al patches at once
%             fitCode=3 %0=fit to mean(Btrans) and mean(gainBtrans)
%                        %1=fit to mean(Btrans) and max(gainBtrans)
%                        %2=fit to max(Btrans) and mean(gainBtrans)
%                        %3=fit to max(Btrans) and max(gainBtrans)
            %[r1,a1,K1,flag1,raR21,r_T1,K_T1,K_T_ratio1,r_T_ratio1]=estSingleSpeciesModel(Btrans,dBtrans,gainBtrans,P,fitCode); %fit growth model to al patches at once
            %fit model to 100% quantiles of biomass and production
            %fitCode=[1 1;0.9 0.9;0.8 0.8];
            fitCode=[.5 .5; .5 .8; .8 .5; .8 .8];
            [r1,a1,z1,K1,flag1,raR21,r_T1,K_T1,K_T_ratio1,r_T_ratio1]=estSingleSpeciesModelmsy(Btrans,dBtrans,gainBtrans,P,fitCode(1,:)); %fit growth model to all patches at once
            %fit model to 90% quantiles of production
            [r2,a2,z2,K2,flag2,raR22,r_T2,K_T2,K_T_ratio2,r_T_ratio2]=estSingleSpeciesModelmsy(Btrans,dBtrans,gainBtrans,P,fitCode(2,:)); %fit growth model to all patches at once
             %fit model to all changes in biomass
            [r3,a3,z3,K3,flag3,raR23,r_T3,K_T3,K_T_ratio3,r_T_ratio3]=estSingleSpeciesModelmsy(Btrans,dBtrans,gainBtrans,P,fitCode(3,:)); %fit growth model to all patches at once
             %fit model to biomass time series
%             [r,a,z,K,flag,raR2,r_T,K_T,K_T_ratio,r_T_ratio]=estSingleSpeciesModel(Btrans,dBtrans,gainBtrans,P,fitCode(4,:)); %fit growth model to all patches and times at once
            [r4,a4,z4,K4,flag4,raR24,r_T4,K_T4,K_T_ratio4,r_T_ratio4]=estSingleSpeciesModelmsy(Btrans,dBtrans,gainBtrans,P,fitCode(4,:)); %fit growth model to all patches at once
%             %fit model to 90% quantiles of biomass and production
%             [r5,a5,z5,K5,flag5,raR25,r_T5,K_T5,K_T_ratio5,r_T_ratio5]=estSingleSpeciesModelmsy(Btrans,dBtrans,gainBtrans,P,fitCode(5,:)); %fit growth model to all patches at once
%              %fit model to 80% quantiles of biomass and production
%             [r6,a6,z6,K6,flag6,raR26,r_T6,K_T6,K_T_ratio6,r_T_ratio6]=estSingleSpeciesModelmsy(Btrans,dBtrans,gainBtrans,P,fitCode(6,:)); %fit growth model to all patches at once
%              %fit model to biomass time series         
            %aEst=nanmean(a.*nanmean(Btrans,3).*raR2./nansum(nanmean(Btrans,3).*raR2));
            %             Bmax=max(Btrans,[],3);
            %             aEst=nanmean(a.*nanmean(Btrans,3)./nansum(nanmean(Btrans,3)));
        end
        
    end
    
    if min(min(min(B)))<log10(eps)
        disp('negative biomass error');
    end
    FoodWebFile=sprintf('Foodweb_%s_%s%s.mat', TimeData, num2str(P.iter), ['_numSpecies' num2str(P.n) '_dT' num2str(P.dT*73000) '_basalSize' num2str(P.s.m0) '_meanD' num2str(HP.sdm(i)) '_pInedible' num2str(P.pInedible) '_fIII']);
    %savemat_Foodweb(FoodWebFile,Z_yrs,B_yrs,gainB_yrs,gainZ_yrs,v_yrs,TE_yrs,PB_yrs,TLik_yrs,TLi_yrs,TLk_yrs,TLall_yrs,Zw_yrs,Bw_yrs,gainBw_yrs,gainZw_yrs,vw_yrs,TEw_yrs,PBw_yrs,TLikw_yrs,TLiw_yrs,TLkw_yrs,TLallw_yrs,BLV_yrs,gainBLV_yrs,BLVw_yrs,gainBLVw_yrs,BLV1_yrs,gainBLV1_yrs,BLV1w_yrs,gainBLV1w_yrs,BLV2_yrs,gainBLV2_yrs,BLV2w_yrs,gainBLV2w_yrs,BLV3_yrs,gainBLV3_yrs,BLV3w_yrs,gainBLV3w_yrs,BLV4_yrs,gainBLV4_yrs,BLV4w_yrs,gainBLV4w_yrs,BLV5_yrs,gainBLV5_yrs,BLV5w_yrs,gainBLV5w_yrs,BLV6_yrs,gainBLV6_yrs,BLV6w_yrs,gainBLV6w_yrs,Btrans,dBtrans,gainBtrans,r,a,z,K,flag,raR2,r_T,K_T,K_T_ratio,r_T_ratio,r1,a1,z1,K1,flag1,raR21,r_T1,K_T1,K_T_ratio1,r_T_ratio1,r2,a2,z2,K2,flag2,raR22,r_T2,K_T2,K_T_ratio2,r_T_ratio2,r3,a3,z3,K3,flag3,raR23,r_T3,K_T3,K_T_ratio3,r_T_ratio3,r4,a4,z4,K4,flag4,raR24,r_T4,K_T4,K_T_ratio4,r_T_ratio4,r5,a5,z5,K5,flag5,raR25,r_T5,K_T5,K_T_ratio5,r_T_ratio5,r6,a6,z6,K6,flag6,raR26,r_T6,K_T6,K_T_ratio6,r_T_ratio6,fitCode,P);
    savemat_Foodweb(FoodWebFile,Z_yrs,B_yrs,gainB_yrs,gainZ_yrs,v_yrs,TE_yrs,PB_yrs,TLik_yrs,TLi_yrs,TLk_yrs,TLall_yrs,Zw_yrs,Bw_yrs,gainBw_yrs,gainZw_yrs,vw_yrs,TEw_yrs,PBw_yrs,TLikw_yrs,TLiw_yrs,TLkw_yrs,TLallw_yrs,BLV1_yrs,gainBLV1_yrs,BLV1w_yrs,gainBLV1w_yrs,BLV2_yrs,gainBLV2_yrs,BLV2w_yrs,gainBLV2w_yrs,BLV3_yrs,gainBLV3_yrs,BLV3w_yrs,gainBLV3w_yrs,BLV4_yrs,gainBLV4_yrs,BLV4w_yrs,gainBLV4w_yrs,Btrans,dBtrans,gainBtrans,r1,a1,z1,K1,flag1,raR21,r_T1,K_T1,K_T_ratio1,r_T_ratio1,r2,a2,z2,K2,flag2,raR22,r_T2,K_T2,K_T_ratio2,r_T_ratio2,r3,a3,z3,K3,flag3,raR23,r_T3,K_T3,K_T_ratio3,r_T_ratio3,r4,a4,z4,K4,flag4,raR24,r_T4,K_T4,K_T_ratio4,r_T_ratio4,fitCode,P);
 
    %plots
    numYrs=20; %number of years in the end of time series to average over
%     plot_demog_spatial_Ball(cat(3,Btrans,B_yrs), P, cat(3,gainBtrans,gainB_yrs),0,numYrs); %latter half of
%     %transcient period + no warming 200 yrs
%     plot_demog_spatial_Ball(cat(3,Btrans,Bw_yrs(:,:,:,end)), P, cat(3,gainBtrans,gainBw_yrs(:,:,:,end)),P.dT(end)*73000,numYrs); %latter half of
% %     %transcient period + warming 200 yrs
% %     plot_demog_spatial_Ball(cat(3,Btrans,Bw_yrs(:,:,:,2)), P, cat(3,gainBtrans,gainBw_yrs(:,:,:,2)),P.dT(2)*73000,numYrs); %latter half of
%     

% plot_demog_spatial_Ball(B_yrs, P, gainB_yrs,0,numYrs);
% plot_demog_spatial_Ball(BLV1_yrs, P, gainBLV1_yrs,0,numYrs);
% plot_demog_spatial_Ball(Btrans, P, gainBtrans,0,numYrs);
%plot_demog_spatial_Ball(cat(3,Btrans,BLV1_yrs), P, cat(3,gainBtrans,gainBLV1_yrs),0,numYrs); %latter half of

end