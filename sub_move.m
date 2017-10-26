function [Bout] = sub_move(Bin)
global P

% FCTS for diffusion
% Bout = ((1-(2.*alpha)).*Bin) + (alpha.*(Bin'*P.Dmat)');
%a = mtimes(Bin',P.Dmat);
a = Bin'*P.Dmat;
%b = bsxfun(@times,P.alpha,a);
b = P.alpha.*a;
c = (1-(2.*P.alpha));
%d = bsxfun(@times,c,Bin');
d = c.*Bin';

Bout = (d + b)';


