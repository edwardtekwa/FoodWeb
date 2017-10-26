function Gij = sub_feedfunc(P);

load ./Data/Data_Barnes
Si = P.S; Sj = P.S; % body sizes
X = [ones(size(P.S')) P.S']; % design matrix for linear regression

% From Barnes et al.
mu = (abs((X*MU.b) - P.S')); % mean difference in prey size (smaller)
sd = (X*VAR.b); % range in prey size difference

%% From Blanchard et al. (alterative source of parameters)
%mu(find(mu)) = 2;
%sd(find(sd)) = 1;

%% Make interaction kernels
for i = 1:length(Si);
 for j = 1:length(Sj);
    if Si(i) - Sj(j) > 0
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

