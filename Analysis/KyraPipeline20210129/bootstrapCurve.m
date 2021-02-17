function [ SEM ] = bootstrapCurve( timesResample,RelevantSamples )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

boots=nans(1,timesResample);
for i=1:timesResample
    boots(i)=nanmean(datasample(RelevantSamples,length(RelevantSamples), 'Replace',true));
end
SEM=std(boots);
end

