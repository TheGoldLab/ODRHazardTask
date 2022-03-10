function [ diff ] = degAngDiff( Ang1, Ang2 )
%UNTITLED4 Summary of this function goes here
%   Given two degrees, finds the difference between the two

diff=rad2deg(angdiff(deg2rad(Ang1),deg2rad(Ang2)));

end

