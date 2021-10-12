function [loc,sca]=gumbel_par(m,s)
% calculate gumbel parameter from mean and standard deviation
sca=sqrt(6)/pi*s;
loc=m-0.5772*sca;
end
