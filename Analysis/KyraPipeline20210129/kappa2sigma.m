function [ sigma ] = kappa2sigma( kappa )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

sigma=180/pi*sqrt(-2*log(besseli(1,kappa)/besseli(0,kappa)));
end

