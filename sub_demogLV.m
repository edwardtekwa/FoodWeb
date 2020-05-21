function [gainB,dB] = sub_demogLV(B,T,r,a,z,P)

% Edward Tekwa Nov 29, 17
% run demographic dynamics based on patch temperatures and estimated
% Lotka-Volterra dynamics

%gainB=B.*skewThEnv(r,T,P.z); %growth
gainB=B.*skewThEnv(r,T,z); %growth
if size(a,1)==2
    dB=gainB+(a(1,:).*B+(nansum(B,2)-B).*a(2,:)).*B; %growth and interactions with self and others
elseif size(a,1)==1
    dB=gainB+(a.*B).*B; %growth and interactions with self only
end

%gainB=(skewThEnv(r,T,z)./-a).*skewThEnv(r,T,z)/4; %overwrite with theoretical maximum production