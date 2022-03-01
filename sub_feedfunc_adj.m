function Gij = sub_feedfunc_adj(P)
%originally by James Watson: based on Barnes et al. 08 dataset
%edited by Edward Tekwa Aug 24, 18: use regression results from Barnes et
%al. 10 instead

%pInedible=portion of possible preys that is inedible (higher means
%specialist predators, with random potential preys knocked out from diet)

%load ./Data/Data_Barnes
Si = P.S; Sj = P.S; % body sizes
X = [ones(size(P.S')) P.S']; % design matrix for linear regression

% From Barnes et al. 08
%mu = (abs((X*MU.b) - P.S')); % mean difference in prey size (smaller)
%sd = (X*VAR.b); % range in prey size difference

% From Barnes et al. 10 (Tekwa)
mu = X*[2.66; 0.24]; % mean difference between predator and prey log10(g)=log10(PPMR)=0.24log10(predator mass)+2.66 (or 2.08, 2.37, 2.95, 3.24)
%sd = (X*[1.36;0]); % range in prey size difference =s.d. of intercept in the log10(PPMR)-log10(predator mass) regression
sd = (X*[0.569;0]);

%% From Blanchard et al. (alterative source of parameters)
%mu(find(mu)) = 2;
%sd(find(sd)) = 1;

%% Make interaction kernels
for i = 1:length(Si); %predators
    for j = 1:length(Sj); %preys
        if Si(i) - Sj(j) > 0 && rand>P.pInedible %if predator is bigger than prey and prey is not selected to be inedible (by probability pInedible)
            a = 1 ./ (sqrt(2.*pi).*sd(i));
            b = exp(-((Si(i)-Sj(j)-mu(i)).^2) ./ (2.*sd(i).^2));
            Gij(i,j) = a .* b;
        else
            Gij(i,j) = 0;
        end
    end
end


% Get rid of zoo feeding
Gij(1,:) = [];

% Normalize to a max of 1
Gij = Gij ./ repmat(max(Gij,[],2),[1 P.n+1]);
%Gij = Gij ./ repmat(sum(Gij,2),[1 P.n+1]); % or normalize to sum to 1

%Gij(Gij<1e-2)=0; %set interactions below threshold to be zero

