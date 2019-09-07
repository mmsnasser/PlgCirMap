function  fz2  = fcaunb(et,etp,f2,z)
% fcaunb.m 
% Nasser, September 6, 2019
%
% This function is a modified version of the MATLAB function "fcaun.m" to
% be used for bounded simply connected domains. 
% 
%
iprec  =  4;
nv2    =  length(et);
w2     =  1/(i*nv2);
vz2    = [real(z) ; imag(z)];       % target
nz2    = length(z);                 % ntarget
a2     = [real(et.') ; imag(et.')]; % source
bf2    = [w2.*f2.*etp].';
[Uf2]  =  zfmm2dpart(iprec,nv2,a2,bf2,0,0,0,nz2,vz2,1,0,0);
b12    = [w2.*etp].';
[U12]  =  zfmm2dpart(iprec,nv2,a2,b12,0,0,0,nz2,vz2,1,0,0);
fz2    = (Uf2.pottarg)./(U12.pottarg);
end