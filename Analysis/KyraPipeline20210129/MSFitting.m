function [ params,err] = MSFitting( LLRs,resp)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% FITTING
resp=resp;
LLRs=LLRs;
model=@fitModel;

%Hrate, noise, lapse)
lb=[0,.1,0];
p0=[.5,1,0];
ub=[1,3,.25];

opts=optimoptions(@fmincon,'Algorithm','interior-point','Display','off');
ms=MultiStart('UseParallel','always');
problem=createOptimProblem('fmincon','objective',...
    model, 'x0', p0, 'lb', lb, 'ub', ub, 'options', opts);
params=run(ms,problem,16);
err=fitModel(params);


function e=fitModel(x)
pc=@(x) PredictedProb(LLRs,x(1),x(2),x(3));
e= -nansum(((1-resp).*log(1-pc(x))+resp.*log(pc(x))));
end

end

