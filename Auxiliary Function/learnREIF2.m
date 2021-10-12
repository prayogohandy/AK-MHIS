function [xstar,maxREIF,idxs,REIF] = learnREIF2(S,ymat,varmat,fp)
% calculate learningn function REIF2
REIF=zeros(size(ymat,1),1);
for i=1:size(ymat,1)
    if mod(i,ceil(size(ymat,1)/100)*10)==0
        fprintf('%f\n',i/size(ymat,1)*100);
    end
    y=ymat(i);sig=sqrt(varmat(i));
    REIF(i,:)=y*(1-2*normcdf(y/sig))+sig*(2-sqrt(2/pi)*exp(-0.5*(y/sig)^2));
end
[maxREIF,idxs]=max(REIF.*fp);
xstar=S(idxs,:);
end