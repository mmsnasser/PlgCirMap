function  fz  = fcauep (et,etp,fep,z,n,finf)
% fcauep.m 
% Nasser, September 6, 2019
%
% This function is a modified version of the MATLAB function "fcau.m" to
% allow us to programe Koebe iterative method for domains with corners.
% 
% 
iprec =  4;
vz    = [real(z) ; imag(z)];       % target
nz    = length(z);                 % ntarget
a     = [real(et.') ; imag(et.')]; % source
tn    = length(et);                % nsource=(m+1)n
bf    = fep.';
[Uf]  = zfmm2dpart(iprec,tn,a,bf,0,0,0,nz,vz,1,0,0);
b1    = [etp].';
[U1]  = zfmm2dpart(iprec,tn,a,b1,0,0,0,nz,vz,1,0,0);
if( nargin == 4 ) 
    fz    = (Uf.pottarg)./(U1.pottarg);
end
if( nargin == 6 ) 
    fz= (finf-(Uf.pottarg)./(n*i))./(1-(U1.pottarg)./(n*i));
end
end