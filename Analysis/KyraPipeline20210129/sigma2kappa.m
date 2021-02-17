function [kappa]=sigma2kappa(sigma)
%Siggma kappa table maker
SKLT=[];
kappa=0:.01:60;
for idx=1:length(kappa)
    SKLT(idx,1)=kappa(idx);
    SKLT(idx,2)=kappa2sigma(kappa(idx));
end

SKLT(:,2)=round(SKLT(:,2),1);

[~,idx]=unique(SKLT(:,2));
SKLT=SKLT(idx,:);
sigma=round(sigma,1);
idx=find(SKLT(:,2) == sigma);
kappa=SKLT(idx,1);
end
