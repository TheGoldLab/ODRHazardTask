function [ PredProb,Belief ] = PredictedProb(LLR,H,ChoiceNoise,lapse )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%Prob of choosing T1, positive belief means pro T1


Belief=PredictedChoices(LLR,H);
PredProb=lapse+(1-2*lapse)./(1+exp(-Belief./ChoiceNoise));
end

