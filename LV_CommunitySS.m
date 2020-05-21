function [SS]=LV_CommunitySS(Bseries,dBseries,Pseries,T,zs,fitCode,params)
%Edward Tekwa Nov 13, 18
%get sum of squares of LV model predictions on Biomass (one species in all patches or(T)emperatures)
%evaluated against data (function of (T)emperature)
%fitCode specifies the: 1-biomass quantile; 2-production quantile, to which
%the model fits (from 0 to 1)

for i=1:length(Bseries,2)
betas=params(1); %maximum r (intrinsic growth rate)
a=params(2); %self-competition
zs=params(3); %optimal temperature for growth (instead of for handling time). Comment out if using handling time thermal optimum


%fit model to a quantile of production and mean of biomass
R=skewThEnv(betas,T,zs); %predicted maximum species intrinsic growth rate at each location
EqB=R/-a; %predicted equilibrium biomass densities at each location
%SSEqB=sum((nanmean(Bseries,3)-EqB).^2); %fit to mean observed biomass
SSEqB=sum((geomean(Bseries,3)-EqB).^2); %fit to geometric mean of observed biomass
%EqGain=skewThEnv(betas,T); %predicted equilibrium productivities at each location
%EqGain=EqB.*skewThEnv(betas,T,zs); %predicted equilibrium species productions at each location
%EqGain=EqB.*skewThEnv(betas,T,zs)/4; %predicted maximum species productions at each location
%SSEqGain=sum((quantile(Pseries,fitCode(2),3)-EqGain).^2); %fit to a quantile of observed growth rate
SSR=sum((quantile(dBseries./Bseries,fitCode(2),3)-R).^2); %fit to a quantile of observed per-capita growth rate
%SS=SSdB/var(dBseries(:))+SSEqB/var(Bseries(:))+SSEqGain/var(Pseries(:));
%SS=SSEqB/var(Bseries(:))+SSEqGain/var(Pseries(:));
SS=SSEqB/var(Bseries(:))+SSR/var(dBseries(:)./Bseries(:));
