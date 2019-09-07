function  fz3  = fcaunu(et,etp,f3,z,nv3,finf)
% fcaunu.m 
% Nasser, September 6, 2019
%
% This function is a modified version of the MATLAB function "fcaun.m" to
% be used for unbounded simply connected domains. 
%
%
iprec  = 4;
w3     = 1/(i*nv3);
vz3    = [real(z) ; imag(z)];       % target
nz3    = length(z);                 % ntarget
a3     = [real(et.') ; imag(et.')]; % source
bf3    = [w3.*f3.*etp].';
[Uf3]  = zfmm2dpart(iprec,nv3,a3,bf3,0,0,0,nz3,vz3,1,0,0);
b13    = [w3.*etp].';
[U13]  = zfmm2dpart(iprec,nv3,a3,b13,0,0,0,nz3,vz3,1,0,0);
fz3= (finf-(Uf3.pottarg))./(1-(U13.pottarg));
end