function  y   =   carg (z)
% carg.m
% Nasser, September 5, 2019
% This function compute the (continuous function) arg(z)
%
%
n       =  length(z);
y       =  unwrap(angle(z));
%
end