function  fz  = fcau(et,etp,f,z,n,finf)
%%
% The function 
%        fz  = fcau(et,etp,f,z,n,finf)
% return the values of the analytic function f computed using the Cauchy
% integral formula at interior vector of points z, where et is the
% parameterization of the boundary, finf is the values of f at infinity 
% for unbounded G, n is the unber of nodes in each boundary component.
% The integral is discretized using the trapezoidal rule. The summations  
% are computed using the FMM.
%
%
%
%
% Please cite this function as:
% M.M.S. Nasser, Fast solution of boundary integral equations with the 
% generalized Neumann kernel, Electronic Transactions on Numerical 
% Analysis,  44 (2015) 189-229.
% 
%
% The MATLAB function fbie depends on the following files:
% 
% zfmm2dpart.m
% fmm2d_r2012a.mexw32
% fmm2d_r2012a.mexw64
% pthreadGC2-w32.dll
% pthreadGC2-w64.dll
% 
% from the Toolbox:
% L. G REENGARD AND Z. G IMBUTAS , FMMLIB2D: A MATLAB toolbox for
% fast multipole method in two dimensions, Version 1.2, 2012.
% http://www.cims.nyu.edu/cmcl/fmm2dlib/fmm2dlib.html
% 
% PLEASE cite the FMMLIB2D toolbox whenever you use the function fcau.m.
%
%
vz    = [real(z) ; imag(z)];       % target
nz    = length(z);                 % ntarget
a     = [real(et.') ; imag(et.')]; % source
tn    = length(et);                % nsource=(m+1)n
iprec = 5;                         %- FMM precision flag
%%
bf    = [f.*etp].';
[Uf]  = zfmm2dpart(iprec,tn,a,bf,0,0,0,nz,vz,1,0,0);
b1    = [etp].';
[U1]  = zfmm2dpart(iprec,tn,a,b1,0,0,0,nz,vz,1,0,0);
if( nargin == 4 ) 
    fz    = (Uf.pottarg)./(U1.pottarg);
end
%%
if( nargin == 6 ) 
    fz= (finf-(Uf.pottarg)./(n*i))./(1-(U1.pottarg)./(n*i));
end
%%
end