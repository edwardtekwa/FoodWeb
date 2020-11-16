function [Prod,dB] = sub_demogLV(B,T,r,a,c,z,Ea,k,smi,Spd) %use NetGainB in output if production excludes intraspecific competition

% Edward Tekwa Nov 29, 17
% run demographic dynamics based on patch temperatures and estimated
% Lotka-Volterra dynamics

%gainB=B.*skewThEnv(r,T,P.z); %growth
gainB=B.*skewThEnv(r,T,z); %max growth


%%%---- Respiration over total population
% (without swimming) rows are species, cols are patches
Met=(0.71.*log(smi) + 18.47 - Ea./(k.*(T+273)))'; %temp dependent
Met=exp(Met);
Met=Met .* (1/7000.*(60.*60.*24)./smi');
% (add swimming) rows are spp, cols are patches
Met = Met .* exp(0.03.*Spd'.*100./60./60./24');
loss_m =  B.* Met';

if size(a,1)==2
    NetGainB=gainB+(a(1,:).*B+(nansum(B,2)-B).*a(2,:)).*B; %growth and interactions with self and others
elseif size(a,1)==1
    NetGainB=gainB+(a.*B).*B; %growth and interactions with self only
end
dB=NetGainB-loss_m;

%Prod=gainB+(a.*B.^2)/c;
%Prod=gainB;
%Prod=(gainB-loss_m)/c;
Prod=gainB-loss_m+(a.*B.^2)./c;
%Prod=loss_m;

%Prod=((skewThEnv(r,T,z)-Met')./-a).*(skewThEnv(r,T,z)-Met')/4; %overwrite with theoretical maximum production