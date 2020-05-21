function [gainB gainZ dB dZ v TE PB] = sub_demog(t,B,Z,T)
global P


% test:
B = B(:,:,1);
Z = Z(:,:,1);
%T = T;


%%---- Define predators and prey
bj = [Z'; B']; % biomass: rows are spp (inc. zooplankton), cols are patches
bj = reshape(bj,[1 P.n+1 P.nx]); % prey biomass: dims are 1,species,patches
bi = reshape(B',[P.n 1 P.nx]); % predator biomass: dims are species,1,patches


%%---- Thermal envelope
% difference between spp optimal temp and local temp: rows patches, cols spp
dt  = bsxfun(@minus, T, P.z);

% squared for thermal envelope
dt2 = dt.^2;

% scaled dt2 for thermal envelope, rows are patches, cols are spp
dt2s= bsxfun(@times, -dt2, 1./P.sigmaz.^2);

% scaled dt for thermal envelope 
dts = bsxfun(@times, P.omegaz./P.sigmaz, dt);

% temp-adjusted search rate, rows patches, cols spp
th    = exp(dt2s).*(1+erf(dts));
th(P.omegaz==-2.7) = th(P.omegaz==-2.7).*0.621525;
th(P.omegaz==0)    = th(P.omegaz==0).*1;
v     = bsxfun(@times, P.V, th); % modified search rate
v     = reshape(v', [P.n 1 P.nx]); % dims: spp, 1, patches


%%---- Type II feeding
%- Prey switching
% Sig = P.Sig; Sig(Sig<0.25) = 0; Sig(Sig>0) = 1; % modulate feeding kernel
% a = bsxfun(@times,bj,Sig).^2; % note exponent (ns)
% b = sum(a,2); c = b(:,ones(P.n+1,1),:); d = a ./ c; % Stock et al. 2008 calculation
% Phi = bsxfun(@times,Sig,d); % note no ms exponent
% Phi(find(isnan(Phi))) = 0; % clean up

Phi = P.Sig; % turn prey switching off

%- encounter rate (numerator of Type II); dims are pred, prey, patch
aa  = bsxfun(@times,v,bj);

%- Factor in prey preference to get numerator of Type II; same dims as aa
bb  = bsxfun(@times,aa,Phi);

%- handling time calculation; dims pred,1 (since all prey),patch
cc  = sum(bsxfun(@times,bb,P.tau),2);
ee  = 1 + cc; % denominator of Type II (ignore interference for now)

% per unit biomass (of predator) predation rate (g/g/day)
cij = bsxfun(@rdivide,bb,ee);
Cij = bsxfun(@times,cij,bi); % get total consumption in g/day: dims are pred, prey, patches



%%---- Zooplankton
loss_z = squeeze(sum(Cij(:,1,:))); % zooplankton that are eaten (g/day/m3)
%gain_z = (P.Zr.*Z) .* (1 - (Z./P.ZK)); % logistic growth in zooplankton (g/day/m3)
gain_z = (P.Zr.*(P.ZK-Z)); % chemostatic growth in zooplankton (g/day/m3), with intrinsic growth rate being the (in and out) flow rate, and carrying capacity being the inflow plankton mass


%%%---- Respiration over total population
% (without swimming) rows are species, cols are patches
%Met = bsxfun(P.m.eq,P.s.mi,T)';
%Met=(exp( (0.71.*log(P.s.mi)) + 18.47 - (0.63./(0.0000862.*(T+273))))...
%                ./7000.*(60.*60.*24)./P.s.mi)';
            
Met=bsxfun(@minus,(0.71.*log(P.s.mi)) + 18.47, (0.63./(0.0000862.*(T+273))))';
Met=exp(Met);
Met=bsxfun(@times,Met,1/7000.*(60.*60.*24)./P.s.mi');


% (add swimming) rows are spp, cols are patches
Met = bsxfun(@times,Met,exp(0.03.*P.Spd'.*100./60./60./24'));

% w/out swimming: rows are spp, cols are patches
loss_m = reshape(bi .* reshape(Met,[P.n 1 P.nx]), [P.n P.nx]);


%%---- Fish consumption and predation
% predation on fish (not on zooplankton) dims are 1, prey, patch
dj = sum(Cij(:,1+1:end,:),1);

% predation mortality (g day-1): rows are prey, cols are patches
loss_f  = reshape(dj,[P.n P.nx]);

% total consumption rate (g day-1): rows are preds, cols are patches
gain_f  = P.lambda .* reshape(sum(Cij,2),[P.n P.nx]);

gainB=gain_f';
gainZ=gain_z;

%%---- Mass balance
dB = gain_f' - (loss_f' + loss_m');
dZ = gain_z - loss_z;

TE = (gain_f - loss_m)./reshape(sum(Cij,2),[P.n P.nx]); % Trophic efficiency
PB = B' ./ (gain_f - loss_m); % biomass/production ratios: days to double

return
