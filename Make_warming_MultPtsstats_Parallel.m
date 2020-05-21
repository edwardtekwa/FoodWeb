%Make_warming_MultPtsstats_Parallel.m
%Edward Tekwa Oct 22, 17
%run food web simulations with parallel warming and no-warming cases
%Nov 3, 17: run multiple temperature change scenarios together and record
%several time points for each simulation, on parallel clusters

%clear all; %close all;
% delete(gcp('nocreate'))
% mycluster=parcluster;
% parpool(min(mycluster.NumWorkers,8)) %run on either the max number of clusters or the limit specified here

numIt=20;    %number of iterations for each parameter combination
SampInt=365;  %record every SampInt pts (days) in time series during transcient period

TimeData=string(datetime);

%timepoints to record:
%TimePts=[100000 118250 136500 154750 173000]; %first point is when temperature change starts; last point is end of simulation (still warming)
%TimePts=[219000:365:292000]; %record every year (600-800 yrs)
RecordYrStart=2000;
TimePts=[RecordYrStart*365+1:365:(RecordYrStart+200)*365+1]; %record every year
%TimePts=[365000:365:438000]; %record every year (1000-1200 yrs)
TransEndYrs=TimePts(1)/2; %number of days at the end of transcient (burn-in) period used to fit (300 yrs)
%single-species model
numPts=length(TimePts); %number of time points
make_traits;
numT=length(tempChange); %number of temperature change scenarios

%parfor
for i = 1:size(TR,2) %go through all parameter combinations (as defined in make_traits.m)
    %% make parameters for a given set of traits
    [P B Z T] = make_parameters(TR,i);
    temps=num2str(P.dT*73000);
    sprintf('Foodweb_%s_%s%s.mat', TimeData, num2str(P.iter), ['_numSpecies' num2str(P.n) '_dT' temps(~isspace(temps)) '_basalSize' num2str(P.s.m0) '_meanD' num2str(HP.sdm(i)) '_stdD' num2str(HP.sdv(i))])
    
    r=0;
    rs=0;
    a=0;
    r1=0;
    a1=0;
    K=0;
    flag=0;
    raR2=0;
    r_T=0;
    K_T_ratio=0;
    r_T_ratio=0;
    
    T1w=[];
    Bw=[];
    Zw=[];
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
    gainBLVw=zeros(P.nx,P.n,numT);
    dBLVw=zeros(P.nx,P.n,numT);
    gainBLV1w=zeros(P.nx,P.n,numT);
    dBLV1w=zeros(P.nx,P.n,numT);
    
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
    BLV_yrs=zeros(P.nx,P.n,numPts);
    gainBLV_yrs=zeros(P.nx,P.n,numPts);
    BLV1_yrs=zeros(P.nx,P.n,numPts);
    gainBLV1_yrs=zeros(P.nx,P.n,numPts);
    
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
    BLVw_yrs=zeros(P.nx,P.n,numPts,numT);
    gainBLVw_yrs=zeros(P.nx,P.n,numPts,numT);
    BLV1w_yrs=zeros(P.nx,P.n,numPts,numT);
    gainBLV1w_yrs=zeros(P.nx,P.n,numPts,numT);
    
    %matrix to record biomasses and changes in species biomasses per patch during transcient period
    dBtrans=zeros(P.nx,P.n,(TimePts(1)-1)/(2*SampInt));
    Btrans=zeros(P.nx,P.n,(TimePts(1)-1)/(2*SampInt));
    gainBtrans=zeros(P.nx,P.n,(TimePts(1)-1)/(2*SampInt));
    
    %% Iterate model forward in time (days at the moment)
    YearStartT=1;
    ttrans=1;
    for t = 1:P.Tend
        if t==TimePts(1) %transition time when no-warming and warming experiments diverge from common states
            %Bw=B; Zw =Z;
            Bw=repmat(B,1,1,numT);
            Zw=repmat(Z,1,1,numT);
            BLV=repmat(B,1,1); %no warming case under estimated Lotka-Volterra dynamics
            BLVw=repmat(B,1,1,numT); %warming case under estimated Lotka-Volterra dynamics
            BLV1=repmat(B,1,1); %no warming case under estimated single-species dynamics
            BLV1w=repmat(B,1,1,numT); %warming case under estimated single-species dynamics
        end
        
        % fix
        B(B<eps) = eps;% eps or 0;
        Z(Z<eps) = eps;% eps;
        if t>=TimePts(1)
            Bw(Bw<eps) = eps;% eps;
            Zw(Zw<eps) = eps;% eps;
            BLV(BLV<eps) = eps;
            BLVw(BLVw<eps)=eps;
            BLV1(BLV1<eps) = eps;
            BLV1w(BLV1w<eps)=eps;
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
            BLV=sub_move(BLV,P); %LV no temp change
            BLV1=sub_move(BLV1,P); %single species model no temp change
            [gainBLV dBLV] = sub_demogLV(BLV,T1,rs,a); %grow/die according to Lotka-Volterra approximation
            [gainBLV1 dBLV1] = sub_demogLV(BLV1,T1,r1,a1);
            for TCase=1:numT %for each temperature change scenario
                Bw(:,:,TCase)        = sub_move(Bw(:,:,TCase),P); % move
                [gainBw(:,:,TCase) gainZw(:,:,TCase) dBw(:,:,TCase) dZw(:,:,TCase) vw(:,:,TCase) TEw(:,:,TCase) PBw(:,:,TCase) TLikw(:,:,TCase) TLiw(:,:,TCase) TLkw(:,:,TCase) TLallw(TCase)] = sub_demog(t,Bw(:,:,TCase),Zw(:,:,TCase),T1w(:,TCase),P); % grow/die
                BLVw(:,:,TCase)        = sub_move(BLVw(:,:,TCase),P); % move
                [gainBLVw(:,:,TCase) dBLVw(:,:,TCase)] = sub_demogLV(BLVw(:,:,TCase),T1w(:,TCase),rs,a); % grow/die
                BLV1w(:,:,TCase)        = sub_move(BLV1w(:,:,TCase),P); % move
                [gainBLV1w(:,:,TCase) dBLV1w(:,:,TCase)] = sub_demogLV(BLV1w(:,:,TCase),T1w(:,TCase),r1,a1); % grow/die

            end
        end
        
        %record transcient time points (every SampInt days):
        if t>TimePts(1)/2 && t<=TimePts(1) &&  floor(t/SampInt)==t/SampInt %last condition is satisfied every SampInt points
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
            BLV_yrs(:,:,tpos)=BLV;
            gainBLV_yrs(:,:,tpos)=gainBLV;
            BLVw_yrs(:,:,tpos,:)=BLVw;
            gainBLVw_yrs(:,:,tpos,:)=gainBLVw;
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
            BLV = BLV + dBLV;
            BLVw = BLVw + dBLVw;
            BLV1 = BLV1 + dBLV1;
            BLV1w = BLV1w + dBLV1w;
        end
        
        %estimate ri(T) and ki(T) for each species at different
        %temperatures, to be used in projection with no species
        %interactions
        if t==TimePts(1)-1
            %estimate per-mass intrinsic growth rate and competition terms
            %[r,a,K,flag,raR2,r_T,rs,K_T_ratio,r_T_ratio]=estSpeciesModel(Btrans,dBtrans,gainBtrans,P); %fit growth model to each patch independently, then estimate intrinsic growth as function of T
            %[r1,a1,K1,flag1,raR21,r_T1,K_T1,K_T_ratio1,r_T_ratio1]=estSingleSpeciesModel(Btrans(end-TransEndYrs:end),dBtrans(end-TransEndYrs:end),gainBtrans(end-TransEndYrs:end),P); %fit growth model to all patches at once
            [r1,a1,K1,flag1,raR21,r_T1,K_T1,K_T_ratio1,r_T_ratio1]=estSingleSpeciesModel(Btrans,dBtrans,gainBtrans,P); %fit growth model to all patches at once
            %aEst=nanmean(a.*nanmean(Btrans,3).*raR2./nansum(nanmean(Btrans,3).*raR2));
            %             Bmax=max(Btrans,[],3);
            %             aEst=nanmean(a.*nanmean(Btrans,3)./nansum(nanmean(Btrans,3)));
        end
        
    end
    
    if min(min(min(B)))<log10(eps)
        disp('negative biomass error');
    end
    FoodWebFile=sprintf('Foodweb_%s_%s%s.mat', TimeData, num2str(P.iter), ['_numSpecies' num2str(P.n) '_dT' num2str(P.dT*73000) '_basalSize' num2str(P.s.m0) '_meanD' num2str(HP.sdm(i)) '_stdD' num2str(HP.sdv(i))]);
    savemat_Foodweb(FoodWebFile,Z_yrs,B_yrs,gainB_yrs,gainZ_yrs,v_yrs,TE_yrs,PB_yrs,TLik_yrs,TLi_yrs,TLk_yrs,TLall_yrs,Zw_yrs,Bw_yrs,gainBw_yrs,gainZw_yrs,vw_yrs,TEw_yrs,PBw_yrs,TLikw_yrs,TLiw_yrs,TLkw_yrs,TLallw_yrs,BLV_yrs,gainBLV_yrs,BLVw_yrs,gainBLVw_yrs,BLV1_yrs,gainBLV1_yrs,BLV1w_yrs,gainBLV1w_yrs,Btrans,dBtrans,gainBtrans,r,a,flag,raR2,rs,K_T_ratio,r_T_ratio,r1,a1,K1,flag1,raR21,r_T1,K_T1,K_T_ratio1,r_T_ratio1,P);
    
    %plots
    %plot_demog_spatial_Btrans(Btrans,Z,prodB,prodZ,P); %latter half of transcient period (no warming)
end