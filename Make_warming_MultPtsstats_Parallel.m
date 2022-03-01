%Make_warming_MultPtsstats_Parallel.m
%Edward Tekwa Oct 22, 17
%run food web simulations with parallel warming and no-warming cases
%Nov 3, 17: run multiple temperature change scenarios together and record
%several time points for each simulation, on parallel clusters
%Dec 19, 21: clean up

delete(gcp('nocreate'))
%% 1. either use the following for parallel runs on the current computer:
%      mycluster=parcluster;
%      parpool(min(mycluster.NumWorkers,10)) %run on either the max number of clusters or the limit specified here
%% 2.or the following for a slurm script to run on a remote cluster (execute the slurm script in terminal to call Make_warming_MultPtsstats_Parallel.m):
   parpool('local', str2num(getenv('SLURM_CPUS_PER_TASK')))
%%
warning('off','all')

numIt=10;    %number of iterations for each parameter combination
SampInt=365;  %record every SampInt pts (days) in time series during transcient period

TimeData=string(datetime); %record time as label for all output files for this simulation set

RecordYrStartRange=[0 0]; %warming starts randomly between the two years [1600 2400]
TransEndYrs=0; %number of years at the end of transcient (burn-in) period used to fit (800)
make_traits; %define parameters to vary
numT=length(tempChange); %number of temperature change scenarios
numYrs=200; %number of warming years (200)
fitCode=[0.5 0.5]; %fit total single species biomass and production to mean food web outcomes during no-warming. Use values other than 0.5 to fit to different quantiles of food web biomass and production. Use 0 to fit model to food web biomass and production at each time point.
numSims=size(TR,2) %total number of simulations
estSimDays=(numSims/str2num(getenv('SLURM_CPUS_PER_TASK')))*8/24 %takes about 8 hrs for 1 simulation per core on Amarel

parfor i = 1:size(TR,2) %go through all parameter combinations (as defined in make_traits.m)
    %% make parameters for a given set of traits
    RecordYrStart=round((RecordYrStartRange(2)-RecordYrStartRange(1))*rand+RecordYrStartRange(1)); %randomize initial time of warming (rounded to day)
    TimePts=[RecordYrStart*365+1:365:(RecordYrStart+numYrs)*365+1]; %record every year, update every day
    numPts=length(TimePts); %number of annual time points to record

    [P B Z] = make_parameters(TR,i); %define parameter values and return initial food web information (P: parameter values, B: species biomass, Z: basal resource biomass - no warming)
    P.Tend=TimePts(end); %set last point to be 200 years after end of random-length transcient period (in days)
    %label for output files:
    FoodWebFile=sprintf('Foodweb_%s_%s%s.mat', TimeData, num2str(P.iter), ['_numSpecies' num2str(P.n) '_dT' num2str(P.dT*73000) '_basalSize' num2str(P.s.m0) '_meanD' num2str(P.sdm) '_pInedible' num2str(P.pInedible) '_fIII'])
    
    %initialize variables:
    %single-species model outputs:
    r1=0;
    a1=0;
    c1=0;
    z1=0;
    K1=0;
    flag1=0;
    raR21=0;
    r_T1=0;
    K_T1=0;
    K_T_ratio1=0;
    r_T_ratio1=0;
    
    %warming (_w) and single-species outputs (_LV). See sub_demog.m for
    %variable definitions
    T1w=[]; %temperature under warming
    Bw=[];
    Zw=[];
    BLV1=[];
    BLV1w=[];
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
    
    %matrices to record time series with no warming:
    Z_yrs=zeros(P.nx,1,numPts);
    B_yrs=zeros(P.nx,P.n,numPts);
    TE_yrs = zeros(P.n, P.nx,numPts);
    PB_yrs = zeros(P.n, P.nx,numPts);
    gainZ_yrs=zeros(P.nx, 1,numPts);
    gainB_yrs=zeros(P.nx, P.n,numPts);
    v_yrs=zeros(P.n,P.nx,numPts);
    TLik_yrs=zeros(P.nx,P.n+1,numPts);
    TLi_yrs=zeros(1,P.n,numPts);
    TLk_yrs=zeros(P.nx,1,numPts);
    TLall_yrs=zeros(numPts,1);
    BLV1_yrs=zeros(P.nx,P.n,numPts);
    gainBLV1_yrs=zeros(P.nx,P.n,numPts);
    
    %matrices to record time series with warming
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
    
    %matrix to record biomasses and changes in species biomasses per patch during transcient period
    dBtrans=zeros(P.nx,P.n,TransEndYrs);
    Btrans=zeros(P.nx,P.n,TransEndYrs);
    gainBtrans=zeros(P.nx,P.n,TransEndYrs);
    %% Iterate model forward in time (days at the moment)
    YearStartT=1;
    ttrans=1;
    for t = 1:P.Tend
        if t==TimePts(1) %transition time when no-warming and warming experiments diverge from common states
            Bw=repmat(B,1,1,numT); %food web species biomass under warming
            Zw=repmat(Z,1,1,numT); %food web basal resource biomass under warming
            BLV1=repmat(B,1,1); %no warming case under estimated single-species dynamics
            BLV1w=repmat(B,1,1,numT); %warming case under estimated single-species dynamics
        end
        
        % fix
        B(B<eps) = 0;
        Z(Z<eps) = 0;
        if t>=TimePts(1)
            Bw(Bw<eps) = 0;
            Zw(Zw<eps) = 0;
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
            [r1,a1,c1,z1,K1,flag1,raR21,r_T1,K_T1,K_T_ratio1,r_T_ratio1]=estSingleSpeciesModel(Btrans,dBtrans,gainBtrans,P,fitCode); %fit growth model to all patches at once
        end
        
    end
    
    if min(min(min(B)))<log10(eps)
        disp('negative biomass error');
    end
    savemat_Foodweb(FoodWebFile,Z_yrs,B_yrs,gainB_yrs,gainZ_yrs,v_yrs,TE_yrs,PB_yrs,TLik_yrs,TLi_yrs,TLk_yrs,TLall_yrs,Zw_yrs,Bw_yrs,gainBw_yrs,gainZw_yrs,vw_yrs,TEw_yrs,PBw_yrs,TLikw_yrs,TLiw_yrs,TLkw_yrs,TLallw_yrs,BLV1_yrs,gainBLV1_yrs,BLV1w_yrs,gainBLV1w_yrs,Btrans,dBtrans,gainBtrans,r1,a1,z1,K1,flag1,raR21,r_T1,K_T1,K_T_ratio1,r_T_ratio1,fitCode,P);
end