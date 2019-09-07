function  [S,Sp]  =  cmb(et,etp,alp,nn,iprec,gmresrestart,gmrestol,gmresmaxit)
% cmb.m
% Nasser, September 5, 2019
% This function implement the numerical method presented in Section 3.2 in
% M.M.S. Nasser, Fast Computation of the Circular Map, Comput. Methods 
% Funct. Theory 15 (2015) 187-223
% to compute the conformal mapping w=f(z) from the bounded simply  
% connected domain onto the interior of the unit circle with the 
% normaliziation f(alpha)=0 and f'(alpha)>0
% 
%
tt         =   (0:2*pi/nn:2*pi-2*pi/nn).';
A          =  et-alp;
Ap         =  etp;
gam  = -log(abs(et-alp));
muh  = -imag(clog(et-alp));
[phi , h ]  =  fbie(et,etp,A,gam,nn,iprec,gmresrestart,gmrestol,gmresmaxit);
S   = phi-muh;
Sp  =  derfft(S-tt)+1;
end