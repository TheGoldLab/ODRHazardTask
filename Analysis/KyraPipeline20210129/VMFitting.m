function [ params,err] = VMFitting( Ev,Mfr,Baseline)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% FITTING
Baseline=Baseline;
Mfr=Mfr;
Ev=Ev;
model=@fitModel;

%Hrate, noise, lapse)
lb=[0,0,3];
p0=[sigma2kappa(15),pi/180*45,10];
ub=[sigma2kappa(60),2*pi,70];

opts=optimoptions(@fmincon,'Algorithm','interior-point','Display','off');
ms=MultiStart('UseParallel','always');
problem=createOptimProblem('fmincon','objective',...
    model, 'x0', p0, 'lb', lb, 'ub', ub, 'options', opts);
params=run(ms,problem,16);
err=fitModel(params);
params=[params,Baseline];


function e=fitModel(x)
yfitM=@(x) x(3)*(exp(x(1)*cos(Ev*pi/180-x(2))))./(2*pi*besseli(0,x(1)))+Baseline;    
e= nansum((Mfr-yfitM(x)).^2);

end

end

