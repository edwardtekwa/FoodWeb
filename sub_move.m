function [Bout] = sub_move(Bin,P)
%global P

% FCTS for diffusion
% Bout = ((1-(2.*alpha)).*Bin) + (alpha.*(Bin'*P.Dmat)');
%a = mtimes(Bin',P.Dmat);
a = Bin'*P.Dmat; %max potential biomass in from neighbouring patches
%b = bsxfun(@times,P.alpha,a);
b = P.alpha.*a; %influx
c = 1-(2.*P.alpha); %fraction of biomass remaining after loss to neighbours (absorbing boundaries)
c1 = 1-P.alpha; %fraction of biomass remaining after loss to neighbours at edge patches (reflecting boundaries)
%d = bsxfun(@times,c,Bin');
d = c.*Bin'; %remaining biomass
%comment out following two lines if using absorbing boundaries
d(:,1)=c1.*Bin(1,:)'; %first patch biomasses lose only 1/2 of other patches due to reflecting boundary
d(:,end)=c1.*Bin(end,:)'; %first patch biomasses lose only 1/2 of other patches due to reflecting boundary


Bout = (d + b)';


