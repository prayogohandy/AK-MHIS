function Pf=impSamp(ysign,fp,hp)
% calculate Pf using importance sampling
Pf=1/length(ysign)*sum(ysign.*fp./hp);
end