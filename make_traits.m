%%%%% Species traits for simulations

% To write out:
% P.n		number of species
% P.nx		number of patches
% P.m.si    body masses (g)
% P.D       diffusivity (m^2/day)
% P.z       thermal optima (deg C)
% P.sigmaz  thermal niche breadth (deg C)
% P.omegaz  thermal niche skewness (deg C)
% P.Binit	initial biomass fish
% P.Zinit	initial biomass zooplankton
% P.sdm	mean D

rng('shuffle') % set random number seed to ensure the same 'random' numbers each run, or shuffle

P.dt=1;
ThermPerfWidth=1.25; %thermal performance width=sqrt(2)*(search performance SD: w_T=0.884)
sdm = [-Inf 0 3 6 7 8]; %[-Inf 0 3 6 7 8]; % mean diffusion coefficient in log10(m^2/day)
sdv = 0; % standard deviation of D among species in log10 space
specialist = [0]; %value of pInedible [0 0.2 0.4 0.6 0.8]
tv = [0]; %standard deviation of optimal temperature around 12.5C: 10~uniform distribution
basalSize =[0.01]; % Body-mass of zooplankton g
tempChange =[3]; %temperature change in degree C
numSpecies=[200]; %number of random initial species
initBDistr=[0]; %initial log body size-log biomass slope

P.nx = 21; % number of patches
minSize = floor(log10(basalSize))+2;
maxSize = 6;
minInitB=log10(eps)+1; %minimum initial biomass in log10 scale
maxInitB=log10(eps)+6; %maximum initial biomass in log10 scale
minOptTemp=0; %minimum optimal temperature
maxOptTemp=34; %maximum optimal temperature

ii=1;
for q=1:length(basalSize)
    for m =1:length(numSpecies)
        for k = 1:length(sdm)
            for j = 1:length(specialist)
                for l = 1:length(initBDistr)
                    for i = 1:numIt % iterations with each hyperparameter combination
                        P.iter=i;
                        P.sdm = sdm(k);
                        P.pInedible =specialist(j);
                        P.n = numSpecies(m); % number of species
                        P.m.si    = 10.^(rand(1,P.n).*(maxSize- minSize)+ minSize); %log uniform random
    
                        % Diffusivity in m^2/day
                        P.D = (10.^(sdm(k)))*ones(1,P.n); %identical diffusion rates for all species
                        %P.D = (10.^(sdm(k))).*P.m.si.^0.46; %diffusion rates increases by power law exponent of 0.46 (derived from Megrey 2007 assuming body length=characteristic length)
                        %P.D = P.D.*10^sdm(k)/mean(P.D); %normalize so that log10(mean(P.D))=sdm(k)
                        
                        % optimal temperature: random in min to max degC
                        P.z  = rand(1,P.n).*(maxOptTemp-minOptTemp)+minOptTemp; %uniform distribution between min and max
                        % thermal niche breadth for species i (deg C, from Urban et al. 2012)
                        P.sigmaz  = ones(1,P.n)*ThermPerfWidth;
                        
                        % skewness of thermal niche for species i (from Urban et al. 2012)
                        P.omegaz  = ones(1,P.n)* -2.7;
                        
                        % initial biomasses rand(1, P.n, 1).*([1:P.n].^0)
                        P.Binit = 10.^(minInitB+(maxInitB-minInitB)*rand(P.nx, P.n, 1)); % random in 10^-5 to 10^0
                        P.ZK = 5; % zooplankton carrying capacity, values from COBALT for g dry weight, times 5g wet fish/1g dry fish
                        P.Zinit = P.ZK;
                        P.s.m0 = basalSize(q); %= 0.05; % Body-mass of zooplankton; g, 5e-7 for phytoplankton, 1e-13 for bacteria
                        P.DT = tempChange; %temperature change in degree C

                        % store
                        TR{ii}    = P;
                        ii = ii + 1;
                    end
                end
            end
        end
    end
end

