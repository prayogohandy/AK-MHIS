function [lam,k]=weibull_par(m,s)
% calculate weibull parameter from mean and standard deviation
f = @(k) s^2/m^2-gamma(1+2./k)./gamma(1+1./k).^2+1;
k0=[0.1,1e10];
k=fzero(f,k0);
lam=m/gamma(1+1/k);
end
