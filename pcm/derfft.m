function y = derfft (x)
% derfft.m
% Nasser, September 6, 2019
%
% Suppose that f(t) is a real 2pi-periodic function on the interval [0,2pi]
% and the 1*n vector x is the values of the function f(t) at the n  
% equidistant points (n must be even)
% t_j=(j-1)*2*pi/n,  j=1,2,...,n.
% The function 
%           y = derfft (x)
% uses fft to find the trigonometric interpolating polynomial that 
% interpolates the function f(t) at the n points t_1,t_2,...,t_n. Then the 
% function derfft computes the values the periodic derivative function 
% f'(t) at the above points t_j.
%
% Example
% n     =  2^7;
% t     = [0:2*pi/n:2*pi-2*pi/n].';
% x     =  sin(t).*exp(-cos(t));
% xpe   =  cos(t).*exp(-cos(t))+sin(t).^2.*exp(-cos(t));
% xp    = derfft (x);
% error = norm(xp-xpe,inf)
%
%
n        =  length(x);
A        =  ifft(x);
%
B        =  zeros(n,1);
for j=1:n/2-1
    B(j+1) = -(1i*j)*A(j+1);
end
for j=n/2+1:n-1
    B(j+1) = -(1i*(j-n))*A(j+1);
end
%
y        =  real(fft(B));
%
end