function [r,a,K,flag,raR2,r_T,K_T,K_T_ratio,r_T_ratio]=estAllSpeciesModel(Btrans,dBtrans,gainBtrans,P) %r, a dimensions: species, patch (temperature)
%Edward Tekwa May 9, 18
%estimate per-mass intrinsic growth rate and competition terms for all
%species at once
%global zs
options=optimoptions(@lsqlin,'OptimalityTolerance',1e-12,'StepTolerance',1e-12,'ConstraintTolerance',1e-12,'display','off');

r=ones(1,size(Btrans,2));
a=ones(1,size(Btrans,2)); %a=NaN(size(Btrans,1),size(Btrans,2)); %first row: self competition, second row: interaction with all other species in patch
raR2=NaN(1,size(Btrans,2)); %model R^2 values for each species
numraPts=NaN(size(Btrans,1),size(Btrans,2));
K=NaN(size(Btrans,1),size(Btrans,2)); %maximum carrying capacity
flag=NaN(1,size(Btrans,2));
%Bmax=max(Btrans,3); %biomass
Bmean=mean(Btrans);
% dBmax=NaN(size(dBtrans,1),size(dBtrans,2)); %change
% dBmean=NaN(size(dBtrans,1),size(dBtrans,2));
%gainBmax=NaN(size(dBtrans,1),size(dBtrans,2)); %productivity
gainBmean=mean(gainBtrans);
numSpecies=size(Btrans,2);
r_T=NaN(size(Btrans,1),size(Btrans,2));
K_T=NaN(size(Btrans,1),size(Btrans,2));
betas=NaN(1,size(Btrans,2));
    
freeparstart=[1e-3 -1e-3]; %r(T),a
freeparmin=[0 -Inf];
freeparmax=[Inf 0];
[lm4,fval4,exitflag,output] = fminsearchbnd(@(params) LV1All_SS(Bmean,gainBmean,P,params),freeparstart,freeparmin,freeparmax); %estimate single-species growth model with skewed thermal envelope for intrinsic growth and temp-independent self-competition
r(:,species)=lm4(1);
a(:,species)=lm4(2);
raR2(species)=1-fval4/(var(Y4)*(length(Y4)-1));
flag(1,species)=exitflag;
%estimate left-skewed normal thermal envelope of species s:
%beta=nlinfit(P.T(nansum(X3(:,1:end-2))>0),lm4(1:end-2),@skewThEnv,1e-5);
r_T(:,species)=skewThEnv(lm4(1),P.T);
%        betas(1,species)=beta;


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