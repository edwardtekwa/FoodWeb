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
% HP.sdm	mean D (hyperparameter)
% HP.sdv	variance of D

rng('shuffle') % set random number seed to ensure the same 'random' numbers each run, or shuffle

%numIt=20;    %number of iterations for each parameter combination
% Experiment #1: Effects of dispersal
% D mean is 10^(___) m^2/day
% D sd is in log10 space
%sdm = [0 2 4 5 7]; % mean D in log10(m^2/day)
%sdv = [0 0.3 0.5 0.7]; % standard deviation of D among species in log10 space
sdm = [-Inf 0 3 6 9]; %[-Inf 0 3 6 9 12]; % mean D in log10(m^2/day) for all
sdv = 0; % standard deviation of D among species in log10 space
specialist = [0 0.5]; %value of pInedible, with extreme specialist->1, and generalist->0
tv = [0]; %standard deviation of optimal temperature around 12.5C: 10~uniform distribution
basalSize =[0.01]; %= 0.05; % Body-mass of zooplankton; g, 5e-7 for phytoplankton, 1e-13 for bacteria
tempChange =[6]; %temperature change in degree C [-2 2 4 6]
%tempChange =[8]; %temperature change in degree C
numSpecies=[200]; %number of random initial species
initBDistr=[0]; %initial log body size-log biomass slope

P.nx = 21; % number of patches
minSize = floor(log10(basalSize))+2; %+3
maxSize = 6;
minInitB=log10(eps)+1; %minimum initial biomass in log10 scale, eg. log(eps) (1)
maxInitB=log10(eps)+6; %maximum initial biomass in log10 scale, eg. log(eps*1e15) (2)
minOptTemp=0; %minimum optimal temperature
maxOptTemp=34; %maximum optimal temperature
%P.Times=TimePts; %time points of warming simulation

ii=1;
%while ii <= length(basalSize); %just iterate through basal size variations for now
for q=1:length(basalSize)
    for m =1:length(numSpecies)
        for k = 1:length(sdm)
            for j = 1:length(specialist)
                for l = 1:length(initBDistr)
                    for i = 1:numIt % iterations with each hyperparameter combination
                        P.iter=i;
                        HP.sdm(ii) = sdm(k);
                        HP.sdv(ii) = sdv;
                        HP.tv(ii) = tv;
                        P.pInedible =specialist(j);
                        P.n = numSpecies(m); % number of species
                        %minSize=floor(log10(basalSize(l)))+1;
                        % body size: random in 10^(-0.4 to 6)
                        %P.m.si    = 10.^(rand(1,P.n).*(6- -0.4)+ -0.4);
                        P.m.si    = 10.^(rand(1,P.n).*(maxSize- minSize)+ minSize); %log uniform random
    
                        % Diffusivity in m^2/day
                        %P.D  = 10.^normrnd(sdm(k), sdv(j), 1, P.n);
                        %P.D  = (P.m.si/(10^maxSize-10^minSize)/2).*(10.^(normrnd(sdm(k), sdv(j), 1, P.n))); %general diffusion function
                        %P.D  = 10.^(log10(P.m.si) .* normrnd(sdm(k), sdv(j), 1, P.n)); %proportional to body size
                        P.D = (10.^(sdm(k)))*ones(1,P.n); %identical diffusion rates for all species
                        
                        % optimal temperature: random in min to max degC
                        P.z  = rand(1,P.n).*(maxOptTemp-minOptTemp)+minOptTemp; %uniform distribution between min and max
                        %P.z = ones(1,P.n).*15;
                        %P.z = (maxSize-minSize-(log10(P.m.si)-minSize))*25/(maxSize-minSize) + normrnd(0, tv(j), 1, P.n); %inversely proportional to body size
                        %P.z = (min(max(ones(1,P.n).*12.5+normrnd(0, tv(j), 1, P.n),minOptTemp),maxOptTemp)); %normal distribution with boundaries at min and max
                        % thermal niche breadth for species i (deg C, from Urban et al. 2012)
                        P.sigmaz  = ones(1,P.n)*5;
                        
                        % skewness of thermal niche for species i (from Urban et al. 2012)
                        P.omegaz  = ones(1,P.n)* -2.7;
                        
                        % initial biomasses rand(1, P.n, 1).*([1:P.n].^0)
                        P.Binit = 10.^(minInitB+(maxInitB-minInitB)*rand(P.nx, P.n, 1)); % random in 10^-5 to 10^0
                        %P.Binit = 10.^(minInitB+(maxInitB-minInitB)*rand(P.nx, P.n, 1)).*([1:P.n].^0); % random in 10^-5 to 10^0
                        %P.ZK = 0.946*5; % zooplankton carrying capacity, values from COBALT for g dry weight, times 5g wet fish/1g dry fish
                        P.ZK = 5;
                        P.Zinit = P.ZK;
                        %P.Zinit = P.ZK*(0.1+(1-0.1)*rand(P.nx, 1, 1)); % fractions of zooplankton carrying capacity from 10 to 100 percent (Ed corrected * from / 1/25/17)
                        P.s.m0 = basalSize(q); %= 0.05; % Body-mass of zooplankton; g, 5e-7 for phytoplankton, 1e-13 for bacteria
                        P.DT = tempChange; %temperature change in degree C
                        %P.DT=0;
                        % store
                        TR{ii}    = P;
                        ii = ii + 1;
                    end
                end
            end
        end
    end
end
% save
save Data/Data_traits.mat TR HP

