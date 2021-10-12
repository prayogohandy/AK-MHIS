function p=jointProb(dist,x)
% calculate joint probability density
PDF=queryPDF(dist,x);
p=prod(PDF')';
end