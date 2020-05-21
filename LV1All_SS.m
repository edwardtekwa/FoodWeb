function [SS]=LV1All_SS(Bmean,gainBmean,P,params)
%get sum of squares of LV model predictions on Biomass (one species in all patches or(T)emperatures)
%evaluated against data (function of (T)emperature)

betas=params(1); %maximum r (intrinsic growth rate)
a=params(2); %self-competition

%[gainB,dB]=sub_demogLV(Bseries,T,betas,a);
%Bguess=Bseries+dB;
EqB=skewThEnv(betas,T,P)/-a; %predicted equilibrium biomass densities at each location
%EqGain=skewThEnv(betas,T); %predicted equilibrium productivities at each location
EqGain=EqB.*skewThEnv(betas,T,P); %predicted equilibrium species productions at each location

%SS=sum((Bseries(2:end)-Bguess(1:end-1)).^2); %compare current observed biomass - previous biomass to predicted dB
%SSdB=sum(sum((dBseries-dB).^2)); %compare current observed change in biomass - predicted dB
%SSEqB=sum((max(sum(Bseries,1))-EqB).^2);
%SSEqGain=sum((max(sum(Pseries,1))-EqGain).^2);
SSEqB=sum((mean(Bseries,3)-EqB).^2);
SSEqGain=sum((mean(Pseries,3)-EqGain).^2);

%SS=SSdB/var(dBseries(:))+SSEqB/var(Bseries(:))+SSEqGain/var(Pseries(:));
SS=SSEqB/var(Bseries(:))+SSEqGain/var(Pseries(:));