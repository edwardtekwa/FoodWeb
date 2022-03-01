function [P B Z] = make_parameters(TR,II)
%%%%% Parameters of size-shifts

%%---- Big parameters
P.iter      = TR{II}.iter;
P.dt        = TR{II}.dt; % time step (days): 1 or 0.1
P.nx        = TR{II}.nx; % number of spatial patches (MUST be an even number)
P.dT        = TR{II}.DT./73000./TR{II}.dt; % rate of temperature increase (deg C per day) over 200 years
P.T         = [linspace(4,24,P.nx)  ]'; %temp gradient
%P.T         = [linspace(15,15,P.nx) ]'; %constant temp

%%---- Zooplankton climatology biomasses (values from COBALT for g dry
%%     weight, times 5g wet fish/1g dry fish
P.ZK        = TR{II}.ZK;
P.Zr        = 0.0075;% g day-1 m-3 zooplankton growth rate

%%---- Initial distribution of traits (size, dispersal variance, thermal envelope)
P.n         = TR{II}.n; % number of species
P.s.m0      = TR{II}.s.m0; % body mass of basal species (g)
P.s.mi      = TR{II}.m.si; % body mass of species (g)
[i j]       = sort(P.s.mi); % sort into ascending order
P.s.mi      = P.s.mi(j);
P.z         = TR{II}.z(j); % Thermal optima (deg C)
P.sigmaz    = TR{II}.sigmaz(j); % Thermal niche breadth (deg C)
P.omegaz    = TR{II}.omegaz(j); % Thermal niche skewness
P.diff      = TR{II}.D(j) .* TR{II}.dt; % log10 diffusion coefficient
P.sdm       = TR{II}.sdm; %mean log10 dispersal diffusion coefficient
P.S         = log10([P.s.m0 P.s.mi]);  % log10 g body sizes of all species (incl. zoo)
P.pInedible = TR{II}.pInedible;


%%---- Consumption parameters
P.lambda    = 0.4 ; % consumption efficiency (gi gj-1) originally: 0.4, 0.2
P.L         = (P.s.mi ./ .012).^ (1./3) ./ 100; % length of individual in m; = (P.s.mi ./ .025).^ (1./3.3) ./ 100;
P.Spd       = (1.*P.s.mi.^0.13)./100.*60.*60.*24; % swim speed from Megrey (mday-1)
P.ga        = 0.26; % fraction of day spent hunting (0.26 or 0.13 or 0.5)
P.V         = P.ga.*pi.*(P.L.^2).*P.Spd./P.s.mi; % max search rate (m3 day-1 g-1)


%%---- movement parameters
% dx in m (explicit FTCS). 
P.Dmat      = diag(ones(1,P.nx-1),1) + diag(ones(1,P.nx-1),-1); %reflective boundaries
%P.Dmat(1,end) = 1; P.Dmat(end,1) = 1; % periodic boundaries
P.alpha     = min(P.diff' .* (1./(((110000/2)*180/P.nx).^2)),1/3); % P.nx grid cells from 90 to 0 degN (90 deg) and 110,000 m per deg N.
%Note: P.alpha maximum is the portion of source biomass that moves to a
%neighbouring patch (=1/3 when there are maximum 2 neighbouring patches)
%eg. Alpha=(10^(logMove))*(1./(((110000/2)*180/P.nx).^2));
%Note: diffusivity that leads to Alpha portion moving to a neighbouring
%patch is
%logMove=log10(Alpha)-log10(1./(((110000/2)*180/P.nx).^2));


%%---- Search rate and handling time parameters
P.T0        = 293.15; % T_0 in K
P.Eh        = -0.71; % Activation energy of handling time
P.Ea        = 0.63; %0.63 for ectotherms or 0.69 for all (0.57, 0.6, 0.63, 0.66, 0.69)
P.k         = 8.6173324e-5; % Boltzmann constant (eV K-1)
%T           = 12.5 + 273; % temperature for handling time in K


%%---- Metabolic costs (day-1) from Brown 2004 as m/Bi=86400*M_hat/(7000s), converting 
P.m.eq = @(s,T) (exp( (0.71.*log(s)) + 18.47 - (P.Ea./(P.k.*(T+273))))...
                ./7000.*(60.*60.*24)./s);             

%%---- vulnerability matrix
P.Sig       = sub_feedfunc_adj(P); % rows are predators, columns are prey (inc. zooplankton)

%%---- Max population growth rate from Brown and/or Ernest 2003 (assuming 12.5 degC)
LogP = (0.761.*log10(P.s.mi./1000)) + 10.85 ...
	- log10(exp(P.Ea./(P.k.*285.5))); % P is kg per ind per yr
P.Pmax = (10.^LogP) ./ P.s.mi .* 1000 ./ 365; %g per g per day

%%---- Handling time (assumming max pop growth from Brown)
%      , in units days gj-1 gi (i is pred, j is prey)
P.tau = @(T) (repmat(reshape((P.lambda ./(((10.^((0.761.*log10(P.s.mi./1000)) + 10.85 - log10(exp(P.Ea./(P.k.*(273+T)))))) ./ P.s.mi .* 1000 ./ 365) + P.m.eq(P.s.mi,T)))',P.n,1,length(T)),[1 P.n+1])); %temp dependent

%%---- Initialize arrays to save final results and interim calculations
% rows are patches, cols are spp, then time
B = zeros(P.nx,P.n); % fish;
Z = zeros(P.nx,1); % zooplankton

% Initial conditions
P.Binit = TR{II}.Binit;
P.Zinit = TR{II}.Zinit;
B(:,:,1) = TR{II}.Binit;
Z(:,1,1) = TR{II}.Zinit;

%%---- clean up
clear bt I J M MU Si Sj VAR X a ans b bv i j mu sd lower hij d sigmaD upper tmax thin
return
