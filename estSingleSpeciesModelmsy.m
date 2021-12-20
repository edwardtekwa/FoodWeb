function [r,a,c,z,K,flag,raR2,r_T,K_T,K_T_ratio,r_T_ratio]=estSingleSpeciesModelmsy(Btrans,dBtrans,gainBtrans,P,fitCode) %r, a dimensions: species, patch (temperature)
%Edward Tekwa Nov 17, 17
%estimate per-mass intrinsic growth rate and competition terms
%global zs
%options=optimoptions(@lsqlin,'OptimalityTolerance',1e-12,'StepTolerance',1e-12,'ConstraintTolerance',1e-12,'display','off');
options = optimset('MaxFunEvals',2000,'MaxIter',2000,'Display','off','TolFun',1e-9,'TolX',1e-9);
%fminoption=optimoptions('display','off');
warning('off','all')

r=NaN(1,size(Btrans,2)); %intrinsic growth rates for all species
a=NaN(1,size(Btrans,2)); %a=NaN(size(Btrans,1),size(Btrans,2)); %first row: self competition, second row: interaction with all other species in patch
c=NaN(1,size(Btrans,2));
z=NaN(1,size(Btrans,2)); %optimal growth temperature for all species
raR2=NaN(1,size(Btrans,2)); %model R^2 values for each species
numraPts=NaN(size(Btrans,1),size(Btrans,2));
K=NaN(size(Btrans,1),size(Btrans,2)); %maximum carrying capacity
flag=NaN(1,size(Btrans,2));
Bmax=NaN(size(Btrans,1),size(Btrans,2)); %biomass
Bmean=NaN(size(Btrans,1),size(Btrans,2));
dBmax=NaN(size(dBtrans,1),size(dBtrans,2)); %change
dBmean=NaN(size(dBtrans,1),size(dBtrans,2));
gainBmax=NaN(size(dBtrans,1),size(dBtrans,2)); %productivity
gainBmean=NaN(size(dBtrans,1),size(dBtrans,2));
r_T=NaN(size(Btrans,1),size(Btrans,2));
K_T=NaN(size(Btrans,1),size(Btrans,2));
betas=NaN(1,size(Btrans,2));

dBperM=dBtrans./Btrans;
for species=1:size(Btrans,2)
    X3=NaN(size(Btrans,1)*size(Btrans,3),size(Btrans,1)+2);
    Y3=NaN(size(Btrans,1)*size(Btrans,3),1); %all dB data for each species
    nPatches=0;
    if sum(Btrans(:,species,end))>eps %check if species is non-extinct at the end of the transcient period
        for patch=1:size(Btrans,1)
            if Btrans(patch,species,end)>eps && max(gainBtrans(patch,species,:))>-Inf
                nPatches=nPatches+1;
                Y=reshape(dBperM(patch,species,:),[],1);
                X1=reshape(Btrans(patch,species,:),[],1); %self biomass within patch (for estimating a_ii)
                X0=reshape(sum(Btrans(patch,[1:species-1,species+1:end],:),2),[],1); %total other biomass within patch (for estimating a_ij)
                Bmax(patch,species)=max(X1);
                Bmean(patch,species)=mean(X1);
                Bmedian(patch,species)=median(X1);
                Bend(patch,species)=X1(end,:);
                dBmax(patch,species)=max(dBtrans(patch,species,:));
                dBmean(patch,species)=mean(dBtrans(patch,species,:));
                gainBmax(patch,species)=max(gainBtrans(patch,species,:));
                gainBmean(patch,species)=mean(gainBtrans(patch,species,:));
                %X2=[ones(size(X1)), X1]; %separate intercepts (r) and biomass-dependence (a) at each location for each species
                %separate intercepts (r) for each location and one biomass-dependence (a)
                X3((patch-1)*length(X1)+1:patch*length(X1),patch)=ones(size(X1));
                X3((patch-1)*length(X1)+1:patch*length(X1),end-1)=X1;
                X3((patch-1)*length(X1)+1:patch*length(X1),end)=X0; %add interaction with others
                Y3((patch-1)*length(X1)+1:patch*length(X1))=Y;
            end
        end
    end
    X4=X3(~isnan(X3(:,end)),:); %delete NaN rows
    X4=X4(:,nansum(X3)>0); %delete NaN columns;
    X4(isnan(X4))=0; %replace NaNs in the design matrix with zero;
    Y4=Y3(~isnan(X3(:,end)),:); %delete NaN rows
    
    if ~isempty(Y4)
        Bseries=Btrans(:,species,:);
        dBseries=dBtrans(:,species,:);
        Pseries=gainBtrans(:,species,:);
        zs=P.z(species); %get current species optimal temperature (for search rate only)
        Eas=P.Ea;
        ks=P.k;
        smis=P.s.mi(species);
        Spds=P.Spd(species);
        freeparstart=[1e-1 -1 1.004]; %r(T),a, and c (aB/c is taken from r as productivity)
        freeparmin=[0 -Inf 1.004]; %best c: 1.004, alternatives: 1.0025, 1.0065
        freeparmax=[Inf 0 1.004];
        [lm4,fval4,exitflag,output] = fminsearchbnd(@(params) LV_SSmsy(Bseries,dBseries,Pseries,P.T,zs,Eas,ks,smis,Spds,fitCode,params),freeparstart,freeparmin,freeparmax,options); %estimate single-species growth model with skewed thermal envelope for intrinsic growth and temp-independent self-competition
        r(:,species)=lm4(1);
        a(:,species)=lm4(2);
        c(:,species)=lm4(3);
        raR2(species)=1-fval4/(var(Y4)*(length(Y4)-1));
        flag(1,species)=exitflag;
        r_T(:,species)=skewThEnv(lm4(1),P.T,zs);
        z(:,species)=zs;
    end
end

K=r./-a(1,:); %max carrying capacities
K_T=r_T./-a(1,:); %carrying capacities based on temperature-dependent r and
sumBmean=nansum(Bmean(:));
sumBmax=nansum(Bmax(:));
sumK=nansum(K(:));
sumK_T=nansum(K_T(:));
K_T_ratio=sumK_T/sumBmean;
r_T_ratio=nansum(K_T(:).*r_T(:)/4)/nansum(gainBmean(:)); %predicted maximum total productivity over observed maximum total productivity