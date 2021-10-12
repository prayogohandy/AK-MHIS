function Y=ballSample(n,d,lb,ub)
% sample n uniform distributed sample inside a d-hyperball with inner
% radius lb and outer radius ub
X=randn(n,d);
lb=lb^d;ub=ub^d;
U=lb+(ub-lb).*rand(n,1);
r=U.^(1/d);
n=vecnorm(X')';
Y=r.*X./n;