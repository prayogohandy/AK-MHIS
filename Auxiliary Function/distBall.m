function hx=distBall(samples,volume)
% calculate the distribution of multi-layered hyperball from cell array of
% samples and volumes
mat=zeros(length(samples));
mat(1,:)=volume/max(volume);
for i=2:length(samples)
    mat(i,1)=size(samples{i},1)/volume(i);
    mat(i,i)=-size(samples{1},1)/volume(1);
end
z=zeros(size(volume));z(1)=1/max(volume);
hx=mat\z;
end
