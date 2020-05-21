function [r,a,z,K,flag,raR2,r_T,K_T,K_T_ratio,r_T_ratio]=estCommunityModel(Btrans,dBtrans,gainBtrans,P,fitCode) %r, a dimensions: species, patch (temperature)
%Edward Tekwa Nov 17, 17
%estimate per-mass intrinsic growth rate and competition terms
%global zs
options=optimoptions(@lsqlin,'OptimalityTolerance',1e-12,'StepTolerance',1e-12,'ConstraintTolerance',1e-12,'display','off');
warning('off','all')

r=NaN(1,size(Btrans,2)); %intrinsic growth rates for all species
a=NaN(1,size(Btrans,2)); %a=NaN(size(Btrans,1),size(Btrans,2)); %first row: self competition, second row: interaction with all other species in patch
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
    %X3=NaN(size(Btrans,1)*size(Btrans,3),size(Btrans,1)+1);
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
                %lm=fitlm(X1,Y, 'y ~ x1');
                %[lm2,resnorm,residual,exitflag,output,lambda]=lsqlin(X2,Y,[-1 -max(X1)],0,[],[],[-Inf;-Inf],[Inf;0],[],options); %r/a<max(Btrans) or x1/x2<=max(Btrans)
                %r(patch,species)=lm.Coefficients.Estimate(1);
                %a(patch,species)=lm.Coefficients.Estimate(2);
                %r(patch,species)=lm2(1);
                %a(patch,species)=lm2(2);
                %flag(patch,species)=exitflag;
                %raR2(patch,species)=lm.Rsquared.Ordinary;
                %numraPts(patch,species)=lm.NumObservations;
            end
        end
    end
    X4=X3(~isnan(X3(:,end)),:); %delete NaN rows
    X4=X4(:,nansum(X3)>0); %delete NaN columns;
    X4(isnan(X4))=0; %replace NaNs in the design matrix with zero;
    Y4=Y3(~isnan(X3(:,end)),:); %delete NaN rows
    %     %construct lower and upper constraints
    %     lower=-Inf*ones(size(X4,2),1);
    %     upper=[Inf*ones(size(X4,2)-1,1);0];
    %construct Bmax constraints
    %nPatches=size(X4,2)-1;
    
    if ~isempty(Y4)
        Bseries=Btrans(:,species,:);
        dBseries=dBtrans(:,species,:);
        Pseries=gainBtrans(:,species,:);
        zs=P.z(species); %get current species optimal temperature (for search rate only)
        freeparstart=[1e-3 -1e-3 12]; %r(T),a, opt growth T
        freeparmin=[0 -Inf -2];
        freeparmax=[Inf 0 40];
        [lm4,fval4,exitflag,output] = fminsearchbnd(@(params) LV_SSmaxRmeanB(Bseries,dBseries,Pseries,P.T,zs,fitCode,params),freeparstart,freeparmin,freeparmax); %estimate single-species growth model with skewed thermal envelope for intrinsic growth and temp-independent self-competition
        r(:,species)=lm4(1);
        a(:,species)=lm4(2);
        z(:,species)=lm4(3);
        raR2(species)=1-fval4/(var(Y4)*(length(Y4)-1));
        flag(1,species)=exitflag;
        %estimate left-skewed normal thermal envelope of species s:
        %beta=nlinfit(P.T(nansum(X3(:,1:end-2))>0),lm4(1:end-2),@skewThEnv,1e-5);
        %r_T(:,species)=skewThEnv(lm4(1),P.T,zs);
        r_T(:,species)=skewThEnv(lm4(1),P.T,lm4(3));
%        betas(1,species)=beta;
    end
end

%prepare initial guesses, minima and maxima of parameters depending on
%whether species has gone extinct
freeparstart=repmat([1e-3 -1e-3 12],1,size(Btrans,2)); %array of parameters (intrinsic growth rates, competition rates, and thermal optima for all species)
freeparmin=repmat([0 -Inf -2],size(Btrans,2),1); %array of parameters (intrinsic growth rates, competition rates, and thermal optima for all species)
freeparmax=repmat([Inf 0 40],size(Btrans,2),1); %array of parameters (intrinsic growth rates, competition rates, and thermal optima for all species)
freeparmin(sum(Btrans(:,:,end))<=eps,:)=0;
freeparmax(sum(Btrans(:,:,end))<=eps,:)=0;
freeparmin(sum(Btrans(:,:,end))<=eps,2)=-Inf;
freeparmax(sum(Btrans(:,:,end))<=eps,2)=-Inf;
freeparmin=reshape(freeparmin,1,[]);
freeparmax=reshape(freeparmax,1,[]);

[lmAll,fvalAll,exitflag,output] = fminsearchbnd(@(params) LV_SSmaxRmeanB(Bseries,dBseries,Pseries,P.T,zs,fitCode,params),freeparstart,freeparmin,freeparmax); %estimate single-species growth model with skewed thermal envelope for intrinsic growth and temp-independent self-competition
r(:,species)=lm4(1);
a(:,species)=lm4(2);
z(:,species)=lm4(3);

%scale r & a from yr^-1 to day^-1
% r=r/365;
% a=a/365;
% r_T=r_T/365;

K=r./-a(1,:); %max carrying capacities
K_T=r_T./-a(1,:); %carrying capacities based on temperature-dependent r and
sumBmean=nansum(Bmean(:));
%sumBmedian=nansum(Bmedian(:));
sumBmax=nansum(Bmax(:));
%sumBend=nansum(Bend(:));
sumK=nansum(K(:));
sumK_T=nansum(K_T(:));
%K_T_ratio=sumK_T/sumBmax %predicted maximum total biomass over observed maximum total biomass
K_T_ratio=sumK_T/sumBmean
%K_Tmedian_ratio=sumK_T/sumBmedian
%K_Tend_ratio=sumK_T/sumBend

r_T_ratio=nansum(K_T(:).*r_T(:))/nansum(gainBmean(:)) %predicted maximum total productivity over observed maximum total productivity

% scrsz = get(0,'ScreenSize');
% figure ('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/5 scrsz(4)]);
% subplot(3,1,1)
% bar(max(r));
% xlabel 'species (ordered by increasing body size)'
% ylabel 'maximum intrinsic growth rate (r)'
% xlim([0,100])
% title([{'Lotka-Volterra parameter estimates from'}; {'food-web simulation'}])
% subplot(3,1,2)
% bar((-a(1,:)).^.1);
% xlabel 'species (ordered by increasing body size)'
% ylabel 'self competition (-a_i_i)^{1/10}'
% xlim([0,100])
% subplot(3,1,3)
% bar(a(2,:));
% xlabel 'species (ordered by increasing body size)'
% ylabel 'competition with others (a_i_o)'
% xlim([0,100])

%temp-independent a