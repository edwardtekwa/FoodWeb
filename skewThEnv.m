function r_T=skewThEnv(max_r,T,z)
%Edward Tekwa Nov 24, 2017
%takes in max_r (the free parameter), Temperature, optimal temp (z) and returns predicted
%intrinsic growth rate (r)

r_T=0.621525*max_r.*(exp(-((T-z-0.4346)/(sqrt(2)*0.8839)).^2)).*(1+erf(-2.7*(T-z-0.4346)/(sqrt(2)*0.8839)));
%0.8839 is the search performance standard deviation
%-2.7 is the skew parameter
%0.621525 is a scaling factor to match r_T at z to be max_r
%-0.4346 is an offset factor so that r_T is maximum at z
