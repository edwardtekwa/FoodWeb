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

rng(10) % set random number seed to ensure the same 'random' numbers each run

numIt=20;    %number of iterations for each parameter combination
% Experiment #1: Effects of dispersal
% D mean is 10^(___) m^2/day
% D sd is in log10 space
%sdm = [0 2 4 5 7]; % mean D in log10(m^2/day)
%sdv = [0 0.3 0.5 0.7]; % standard deviation of D among species in log10 space
sdm = [-Inf,0,1,10]; % mean D in log10(m^2/day) for all
sdv = [0]; % standard deviation of D among species in log10 space
tv = [0]; %standard deviation of optimal temperature around 12.5C: 10~uniform distribution
basalSize =[0.01]; %= 0.05; % Body-mass of zooplankton; g, 5e-7 for phytoplankton, 1e-13 for bacteria

P.n = 100; % number of species
P.nx = 11; % number of patches
minSize = floor(log10(basalSize))+1;
maxSize = 6;
minInitB=log10(eps)+1; %minimum initial biomass in log10 scale, eg. log(eps)
maxInitB=log10(eps)+2; %maximum initial biomass in log10 scale, eg. log(eps*1e15)
minOptTemp=0; %minimum optimal temperature
maxOptTemp=25; %maximum optimal temperature

ii=1;
while ii <= length(basalSize); %just iterate through basal size variations for now
    for k = 1:length(sdm) % for each hyperparameter value
        for j = 1:length(sdv)
            for l = 1:length(basalSize)
                for i = 1:numIt % iterations with each hyperparameter combination
                    HP.sdm(ii) = sdm(k);
                    HP.sdv(ii) = sdv(j);
                    HP.tv(ii) = tv;
                    %minSize=floor(log10(basalSize(l)))+1;
                    % body size: random in 10^(-0.4 to 6)
                    %P.m.si    = 10.^(rand(1,P.n).*(6- -0.4)+ -0.4);
                    P.m.si    = 10.^(rand(1,P.n).*(maxSize- minSize)+ minSize);
                    
                    % Diffusivity in m^2/day
                    P.D  = 10.^normrnd(sdm(k), sdv(j), 1, P.n);
                    %P.D  = 10.^(log10(P.m.si) .* normrnd(sdm(k), sdv(j), 1, P.n)); %proportional to body size
                    
                    % optimal temperature: random in min to max degC
                    P.z  = rand(1,P.n).*(maxOptTemp-minOptTemp)+minOptTemp; %uniform distribution between min and max
                    %P.z = ones(1,P.n).*15;
                    %P.z = (maxSize-minSize-(log10(P.m.si)-minSize))*25/(maxSize-minSize) + normrnd(0, tv(j), 1, P.n); %inversely proportional to body size
                    %P.z = (min(max(ones(1,P.n).*12.5+normrnd(0, tv(j), 1, P.n),minOptTemp),maxOptTemp)); %normal distribution with boundaries at min and max
                    % thermal niche breadth for species i (deg C, from Urban et al. 2012)
                    P.sigmaz  = ones(1,P.n)*5;
                    
                    % skewness of thermal niche for species i (from Urban et al. 2012)
                    P.omegaz  = ones(1,P.n)* -2.7;
                    
                    % initial biomasses
                    P.Binit = 10.^(minInitB+(maxInitB-minInitB)*rand(P.nx, P.n, 1)); % random in 10^-5 to 10^0
                    P.ZK = 0.946*5; % zooplankton carrying capacity, values from COBALT for g dry weight, times 5g wet fish/1g dry fish
                    P.Zinit = P.ZK*(0.1+(1-0.1)*rand(P.nx, 1, 1)); % fractions of zooplankton carrying capacity from 10 to 100 percent (Ed corrected * from / 1/25/17)
                    P.s.m0 = basalSize(l); %= 0.05; % Body-mass of zooplankton; g, 5e-7 for phytoplankton, 1e-13 for bacteria
                    
                    % store
                    TR{ii}    = P;
                    ii = ii + 1;
                end
            end
        end
    end
end
% save
save Data/Data_traits.mat TR HP

