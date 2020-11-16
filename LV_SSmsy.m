function [SS]=LV_SSmsy(Bseries,dBseries,Pseries,T,zs,Eas,ks,smis,Spds,fitCode,params)
%Edward Tekwa Aug 24, 18
%get sum of squares of LV model predictions on Biomass (one species in all patches or(T)emperatures)
%evaluated against data (function of (T)emperature)
%fitCode specifies the: 1-biomass quantile; 2-production quantile, to which
%the model fits (from 0 to 1)

r=params(1); %maximum r (intrinsic growth rate)
a=params(2); %self-competition
c=params(3);

%zs=params(4); %optimal temperature for growth (instead of for handling time). Comment out if using handling time thermal optimum
%%%---- Respiration over total population
% (without swimming) rows are species, cols are patches
Met=(0.71.*log(smis) + 18.47 - Eas./(ks.*(T+273)))'; %temp dependent
Met=exp(Met);
Met=Met .* (1/7000.*(60.*60.*24)./smis');
% (add swimming) rows are spp, cols are patches
Met = Met .* exp(0.03.*Spds'.*100./60./60./24');
%loss_m =  Bseries.* Met';

if fitCode(1)==0 %fit model to all observed biomass time points
    [gainB,dB]=sub_demogLV(Bseries,T,r,a,c,zs,Eas,ks,smis,Spds); %B,T,r,a,z,P
    %Bguess=Bseries+dB;
    SS=sum(sum((dBseries-dB).^2)); %compare current observed biomass - previous biomass to predicted dB for each time and location
    %SS=sum((Bseries(2:end)-Bguess(1:end-1)).^2); %compare current observed biomass - previous biomass to predicted dB
    %SSdB=sum(sum((dBseries-dB).^2)); %compare current observed change in biomass - predicted dB

else %fit model to a quantile of biomass and production   
    EqB=(skewThEnv(r,T,zs)-Met')/-a; %predicted equilibrium biomass densities at each location
    if fitCode(1)==0.5 %use means
        SSEqB=sum((nanmean(Bseries,3)-EqB).^2); %fit to mean of observed biomass
    else
        SSEqB=sum((quantile(Bseries,fitCode(1),3)-EqB).^2); %fit to a quantile of observed biomass
    end
    %EqGain=EqB.*skewThEnv(r,T,zs); %predicted equilibrium species productions at each location
    %EqGain=EqB.*(skewThEnv(r,T,zs)-Met'); %predicted equilibrium species productions at each location
    %EqGain=EqB.*(skewThEnv(r,T,zs)-Met')/4; %predicted maximum species productions at each location based on maximum sustainable yield
    %EqGain=EqB.*(skewThEnv(r,T,zs)+(a.*EqB)); %predicted maximum species productions at each location based on production at equilibrium biomass (carrying capacity)
    %EqGain=EqB.*(skewThEnv(r,T,zs)+(a.*EqB)/c); %predicted maximum species productions at each location based on production at equilibrium biomass (carrying capacity)
    EqGain=EqB.*(skewThEnv(r,T,zs)-Met'+(a.*EqB)./c); %predicted maximum species productions at each location based on production at equilibrium biomass (carrying capacity)
    %EqGain=EqB.*Met'; %=EqB.*(skewThEnv(r,T,zs)+(a.*EqB)) at equilibrum
    %EqGain=EqB.*(skewThEnv(r,T,zs)-Met')/c;
    
    if fitCode(2)==0.5 %use means
        SSEqGain=sum((nanmean(Pseries,3)-EqGain).^2); %fit to mean of observed production 
    else
        SSEqGain=sum((quantile(Pseries,fitCode(2),3)-EqGain).^2); %fit to a quantile of observed production
    end
    %SS=SSdB/var(dBseries(:))+SSEqB/var(Bseries(:))+SSEqGain/var(Pseries(:));
    %SS=SSEqB/var(Bseries(:))+SSEqGain/var(Pseries(:));
%     if nanmean(Bseries(:))==0
%         ratio_prod_biom=1;
%     elseif isnan(nanmean(Bseries(:)))
%         ratio_prod_biom=1;
%     else
%         ratio_prod_biom=nanmean(Pseries(:))/nanmean(Bseries(:));
%     end
    %SS=SSEqB/abs(Bseries(:))+SSEqGain/abs(nanmean(Pseries(:)));
    SS=SSEqB+1000*SSEqGain;
    %SS=SSEqB*SSEqGain;
end