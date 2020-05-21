function r_T=skewThEnv(max_r,T,z)
%Edward Tekwa Nov 24, 2017
%takes in max_r (the free parameter), Temperature, optimal temp (zs) and returns predicted
%intrinsic growth rate (r)
%global zs

r_T=0.621525*max_r.*(exp(-((T-z-1.738211)/5).^2)).*(1+erf(-2.7*(T-z-1.738211)/5));