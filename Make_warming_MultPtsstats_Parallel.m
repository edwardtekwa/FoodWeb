%Make_warming_MultPtsstats_Parallel.m
%Edward Tekwa Oct 22, 17
%run food web simulations with parallel warming and no-warming cases
%Nov 3, 17: run multiple temperature change scenarios together and record
%several time points for each simulation, on parallel clusters
%Dec 19, 21: clean up

delete(gcp('nocreate'))
%% either use the following for parallel runs on the current computer:
% mycluster=parcluster;
% parpool(min(mycluster.NumWorkers,10)) %run on either the max number of clusters or the limit specified here
%%or the following for a slurm script to run on a remote cluster (execute the slurm script in terminal to call this file):
parpool('local', str2num(getenv('SLURM_CPUS_PER_TASK')))
%%
warning('off','all')

numIt=10;    %number of iterations for each parameter combination: variations: pIned=0.1, alpha_R=2.08, Ea=0.69, Fh=0.13, lambda=0.2
SampInt=365;  %record every SampInt pts (days) in time series during transcient period

TimeData=string(datetime);


RecordYrStartRange=[1600 2400] %warming starts randomly between the two years 1600-2400
TransEndYrs=800; %number of years at the end of transcient (burn-in) period used to fit 800
make_traits;
numT=length(tempChange); %number of temperature change scenarios

numSims=numIt*length(specialist)*length(sdm)
estSimDays=(numSims/str2num(getenv('SLURM_CPUS_PER_TASK')))*12/24 %takes about 12 hrs for 1 simulation per core on Amarel

parfor i = 1:size(TR,2) %go through all parameter combinations (as defined in make_traits.m)
    %% make parameters for a given set of traits
    RecordYrStart=round((RecordYrStartRange(2)-RecordYrStartRange(1))*rand+RecordYrStartRange(1)); %randomize initial time of warming (rounded to day)
    TimePts=[RecordYrStart*365+1:365:(RecordYrStart+200)*365+1]; %record every year
    numPts=length(TimePts); %number of time points

    [P B Z T] = make_parameters(TR,i);
    P.Tend=TimePts(end); %set last point to be 200 years after end of random-length transcient period (in days)
    temps=num2str(P.dT*73000);
    sprintf('Foodweb_%s_%s%s.mat', TimeData, num2str(P.iter), ['_numSpecies' num2str(P.n) '_dT' temps(~isspace(temps)) '_basalSize' num2str(P.s.m0) '_meanD' num2str(HP.sdm(i)) '_stdD' num2str(HP.sdv(i))])

    r1=0;
    a1=0;
    c1=0;
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
    K=0;
    flag=0;
    raR2=0;
    r_T=0;
    K_T_ratio=0;
    r_T_ratio=0;
    
    T1w=[];
    Bw=[];
    Zw=[];
    BLV1=[];
    BLV1w=[];
    BLV2=[];
    BLV2w=[];
    BLV3=[];
    BLV3w=[];
    BLV4=[];
    BLV4w=[];
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
    gainBLV1w=zeros(P.nx,P.n,numT);
    dBLV1w=zeros(P.nx,P.n,numT);
    gainBLV2w=zeros(P.nx,P.n,numT);
    dBLV2w=zeros(P.nx,P.n,numT);
    gainBLV3w=zeros(P.nx,P.n,numT);
    dBLV3w=zeros(P.nx,P.n,numT);
    gainBLV4w=zeros(P.nx,P.n,numT);
    dBLV4w=zeros(P.nx,P.n,numT);
    
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
    BLV1_yrs=zeros(P.nx,P.n,numPts);
    gainBLV1_yrs=zeros(P.nx,P.n,numPts);
    BLV2_yrs=zeros(P.nx,P.n,numPts);
    gainBLV2_yrs=zeros(P.nx,P.n,numPts);
    BLV3_yrs=zeros(P.nx,P.n,numPts);
    gainBLV3_yrs=zeros(P.nx,P.n,numPts);
    BLV4_yrs=zeros(P.nx,P.n,numPts);
    gainBLV4_yrs=zeros(P.nx,P.n,numPts);
    
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
    BLV1w_yrs=zeros(P.nx,P.n,numPts,numT);
    gainBLV1w_yrs=zeros(P.nx,P.n,numPts,numT);
    BLV2w_yrs=zeros(P.nx,P.n,numPts,numT);
    gainBLV2w_yrs=zeros(P.nx,P.n,numPts,numT);
    BLV3w_yrs=zeros(P.nx,P.n,numPts,numT);
    gainBLV3w_yrs=zeros(P.nx,P.n,numPts,numT);
    BLV4w_yrs=zeros(P.nx,P.n,numPts,numT);
    gainBLV4w_yrs=zeros(P.nx,P.n,numPts,numT);
    
    %matrix to record biomasses and changes in species biomasses per patch during transcient period
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
            BLV1=repmat(B,1,1); %no warming case under estimated single-species dynamics
            BLV1w=repmat(B,1,1,numT); %warming case under estimated single-species dynamics
            BLV2=repmat(B,1,1); %no warming case under estimated single-species dynamics
            BLV2w=repmat(B,1,1,numT); %warming case under estimated single-species dynamics
            BLV3=repmat(B,1,1); %no warming case under estimated single-species dynamics
            BLV3w=repmat(B,1,1,numT); %warming case under estimated single-species dynamics
            BLV4=repmat(B,1,1); %no warming case under estimated single-species dynamics
            BLV4w=repmat(B,1,1,numT); %warming case under estimated single-species dynamics
        end
        
        % fix
        B(B<eps) = 0;% eps or 0;
        Z(Z<eps) = 0;% eps;
        if t>=TimePts(1)
            Bw(Bw<eps) = 0;% eps;
            Zw(Zw<eps) = 0;% eps;
            BLV1(BLV1<eps) = 0;
            BLV1w(BLV1w<eps)= 0;
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
            [gainBLV1 dBLV1] = sub_demogLV(BLV1,T1,r1,a1,c1,z1,P.Ea,P.k,P.s.mi,P.Spd);
            for TCase=1:numT %for each temperature change scenario
                Bw(:,:,TCase)        = sub_move(Bw(:,:,TCase),P); % move
                [gainBw(:,:,TCase) gainZw(:,:,TCase) dBw(:,:,TCase) dZw(:,:,TCase) vw(:,:,TCase) TEw(:,:,TCase) PBw(:,:,TCase) TLikw(:,:,TCase) TLiw(:,:,TCase) TLkw(:,:,TCase) TLallw(TCase)] = sub_demog(t,Bw(:,:,TCase),Zw(:,:,TCase),T1w(:,TCase),P); % grow/die
                BLV1w(:,:,TCase)        = sub_move(BLV1w(:,:,TCase),P); % move
                [gainBLV1w(:,:,TCase) dBLV1w(:,:,TCase)] = sub_demogLV(BLV1w(:,:,TCase),T1w(:,TCase),r1,a1,c1,z1,P.Ea,P.k,P.s.mi,P.Spd); % grow/die
            end
        end
        
        %record transcient time points (every SampInt days):
        if t>(TimePts(1)-TransEndYrs*365) && t<=TimePts(1) &&  floor(t/SampInt)==t/SampInt %last condition is satisfied every SampInt points
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
            BLV1_yrs(:,:,tpos)=BLV1;
            gainBLV1_yrs(:,:,tpos)=gainBLV1;
            BLV1w_yrs(:,:,tpos,:)=BLV1w;
            gainBLV1w_yrs(:,:,tpos,:)=gainBLV1w;
        end
        
        % update biomasses
        Z = Z + (dZ * P.dt);
        B = B + (dB * P.dt);
        if t>=TimePts(1)
            Zw = Zw + (dZw * P.dt);
            Bw = Bw + (dBw * P.dt);
            BLV1 = BLV1 + dBLV1;
            BLV1w = BLV1w + dBLV1w;
        end
        
        %estimate ri(T) and ki(T) for each species at different
        %temperatures, to be used in projection with no species
        %interactions
        if t==TimePts(1)-1
            fitCode=[0.5 0.5]; %fit single species model to mean biomass and production in no-warming period
            [r1,a1,c1,z1,K1,flag1,raR21,r_T1,K_T1,K_T_ratio1,r_T_ratio1]=estSingleSpeciesModelmsy(Btrans,dBtrans,gainBtrans,P,fitCode(1,:)); %fit growth model to all patches at once
        end
        
    end
    
    if min(min(min(B)))<log10(eps)
        disp('negative biomass error');
    end
    FoodWebFile=sprintf('Foodweb_%s_%s%s.mat', TimeData, num2str(P.iter), ['_numSpecies' num2str(P.n) '_dT' num2str(P.dT*73000) '_basalSize' num2str(P.s.m0) '_meanD' num2str(HP.sdm(i)) '_pInedible' num2str(P.pInedible) '_fIII']); %pIned=0.1, alpha_R=2.08, Ea=0.69, Fh=0.13, lambda=0.2
    savemat_Foodweb(FoodWebFile,Z_yrs,B_yrs,gainB_yrs,gainZ_yrs,v_yrs,TE_yrs,PB_yrs,TLik_yrs,TLi_yrs,TLk_yrs,TLall_yrs,Zw_yrs,Bw_yrs,gainBw_yrs,gainZw_yrs,vw_yrs,TEw_yrs,PBw_yrs,TLikw_yrs,TLiw_yrs,TLkw_yrs,TLallw_yrs,BLV1_yrs,gainBLV1_yrs,BLV1w_yrs,gainBLV1w_yrs,Btrans,dBtrans,gainBtrans,r1,a1,z1,K1,flag1,raR21,r_T1,K_T1,K_T_ratio1,r_T_ratio1,fitCode,P);
end