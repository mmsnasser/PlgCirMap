function  fz  = fcaun(et,etp,f,z,nv,finf)
% fcaun.m 
% Nasser, September 6, 2019
%
% This function is a modified version of the MATLAB function "fcau.m" to
% allow us to use different numbers of nodes in the boundary components
% of the domain. Recall that nv is a vector of length m where nv(k) is the 
% number of nodes on the boundary component \Gamma_k, k=1,2,...,m. 
% 
% 
% if( nargin == 4 ) 
%     nv    = length(et); % For simply connected domains
% end
m        =  length(nv);
for k=1:m
    w(1+sum(nv(1:k-1)):sum(nv(1:k)),1)=1/(i*nv(k));
end
vz    = [real(z) ; imag(z)];       % target
nz    = length(z);                 % ntarget
a     = [real(et.') ; imag(et.')]; % source
tn    = length(et);                % nsource=el
iprec = 4;                         %- FMM precision flag
bf    = [w.*f.*etp].';
[Uf]  = zfmm2dpart(iprec,tn,a,bf,0,0,0,nz,vz,1,0,0);
b1    = [w.*etp].';
[U1]  = zfmm2dpart(iprec,tn,a,b1,0,0,0,nz,vz,1,0,0);
if( nargin == 4 | 5 ) 
    fz    = (Uf.pottarg)./(U1.pottarg);
end
if( nargin == 6 ) 
    fz= (finf-(Uf.pottarg))./(1-(U1.pottarg));
end
end