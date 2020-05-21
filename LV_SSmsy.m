function [SS]=LV_SSmsy(Bseries,dBseries,Pseries,T,zs,fitCode,params)
%Edward Tekwa Aug 24, 18
%get sum of squares of LV model predictions on Biomass (one species in all patches or(T)emperatures)
%evaluated against data (function of (T)emperature)
%fitCode specifies the: 1-biomass quantile; 2-production quantile, to which
%the model fits (from 0 to 1)

betas=params(1); %maximum r (intrinsic growth rate)
a=params(2); %self-competition
%zs=params(3); %optimal temperature for growth (instead of for handling time). Comment out if using handling time thermal optimum

if fitCode(1)==0 %fit model to all observed biomass time points
    [gainB,dB]=sub_demogLV(Bseries,T,betas,a,zs,[]); %B,T,r,a,z,P
    %Bguess=Bseries+dB;
    SS=sum(sum((dBseries-dB).^2)); %compare current observed biomass - previous biomass to predicted dB for each time and location
    %SS=sum((Bseries(2:end)-Bguess(1:end-1)).^2); %compare current observed biomass - previous biomass to predicted dB
    %SSdB=sum(sum((dBseries-dB).^2)); %compare current observed change in biomass - predicted dB

else %fit model to a quantile of biomass and production   
    EqB=skewThEnv(betas,T,zs)/-a; %predicted equilibrium biomass densities at each location
    %if fitCode(1)==0.5 %use means
    %    SSEqB=sum((nanmean(Bseries,3)-EqB).^2); %fit to mean of observed biomass
    %else
        SSEqB=sum((quantile(Bseries,fitCode(1),3)-EqB).^2); %fit to a quantile of observed biomass
    %end
    %EqGain=skewThEnv(betas,T); %predicted equilibrium productivities at each location
    %EqGain=EqB.*skewThEnv(betas,T,zs); %predicted equilibrium species productions at each location
    EqGain=EqB.*skewThEnv(betas,T,zs)/4; %predicted maximum species productions at each location
    %if fitCode(2)==0.5 %use means
    %    SSEqGain=sum((nanmean(Pseries,3)-EqGain).^2); %fit to mean of observed production 
    %else
        SSEqGain=sum((quantile(Pseries,fitCode(2),3)-EqGain).^2); %fit to a quantile of observed production
    %end
    %SS=SSdB/var(dBseries(:))+SSEqB/var(Bseries(:))+SSEqGain/var(Pseries(:));
    SS=SSEqB/var(Bseries(:))+SSEqGain/var(Pseries(:));
end