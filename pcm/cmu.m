function  [S,Sp,c]  =  cmu(et,etp,alp,nn,iprec,gmresrestart,gmrestol,gmresmaxit)
% cmu.m
% Nasser, September 5, 2019
% This function implement the numerical method presented in Section 3.3 in
% M.M.S. Nasser, Fast Computation of the Circular Map, Comput. Methods 
% Funct. Theory 15 (2015) 187-223
% to compute the conformal mapping w=f(z) from the unbounded simply  
% connected domain onto the exterior of the unit circle with the 
% normaliziation f(z)=z+O(1/z) near inf
%
%
tt        =   (0:2*pi/nn:2*pi-2*pi/nn).';
A          =  ones(nn,1);
Ap         =  zeros(nn,1);
gam  =  log(abs(et-alp));
muh  =  imag(clog(et-alp));
[phi , h ]  =  fbie(et,etp,A,gam,nn,iprec,gmresrestart,gmrestol,gmresmaxit);
S   = phi-muh;
ho  = sum(h)/nn;
Sp  =  derfft(S-tt)+1;
c         =  exp(+ho);
end