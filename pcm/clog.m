function  y   =   clog (z)
% clog.m
% Nasser, September 5, 2019
% This function compute the a continuous branch function of log(z)
%
%
%%
y     =     log(abs(z))+i.*carg(z);
%%