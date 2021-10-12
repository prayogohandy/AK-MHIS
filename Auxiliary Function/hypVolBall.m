function hypvol=hypVolBall(dim,inrad,outrad)
% calculate hypervolume of d-hyperball with inner radius inrad and outer
% radius outrad
hypvol=pi^(dim/2)/gamma(dim/2+1)*(outrad^dim-inrad^dim);
end