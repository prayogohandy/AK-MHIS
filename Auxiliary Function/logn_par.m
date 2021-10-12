function [mu,sig]=logn_par(m,s)
% calculate lognormal parameter from mean and standard deviation
mu=log(m^2/sqrt(s^2+m^2));
sig=sqrt(log(s^2/m^2+1));
end
