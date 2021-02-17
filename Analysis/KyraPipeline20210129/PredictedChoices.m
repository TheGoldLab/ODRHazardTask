function [ Beliefs ] = PredictedChoices(LLRs,H)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% negative means T1 favored, pos means T2 favored)

Beliefs=nans(size(LLRs));
Beliefs(1)=LLRs(1);

for i=2:length(Beliefs)
 Beliefs(i)=LLRs(i)+Beliefs(i-1)+log((1-H)/H+exp(-Beliefs(i-1)))-log((1-H)/H+exp(Beliefs(i-1)));
%disp(PredictChoices(i));
end
% Beliefs(Beliefs==0)=.00000001;

end

