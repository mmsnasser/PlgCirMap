function   f = plgcirmap(ver,alpha,preimg1)
% plgcirmap.m (polygons domains to circular domains map)
% Nasser, September 3, 2019
%
% This MATLAB function compute the conformal mapping w=f(z) from a polygon
% domain G onto a circular domin D and its invers z=f^-1(w). The domain G
% is multiply connected of connectivity m with the boundary
% \Gamma=\Gamma_1 U ... U \Gamma_m. If G is bounded, then the external
% boundary component is \Gamma_m.
% When m=1, the domain G is simply connected domain.
% 
% Input:
% ver (ver is a Cell Array where ver{k} contains the vertices of the 
%      polygon k, k=1,2,...,m).
% alpha (alpha is an auxiliary point in the domain G if G is bounded
%        and alpha=inf if G is unbounded) 
% preimg1 (optional) (preimg1 is the last vertix of the last polygon)
%          If preimg1 is given, then the conformal mapping is normalized
%          by:  f(alpha)=0 and f(ver{end}(end))=1
%          If preimg1 is not given, then the conformal mapping is
%          normalized by:  f(alpha)=0 and f'(alpha)>0
% 
% Output: 
%     f, an object with:
% 
% f.ver = ver (ver is as above)
% f.alpha = alpha (alpha is as above)
% f.nv = nv (nv is a vector of length m where nv(k) is the number of nodes
%            on the boundary component \Gamma_k, k=1,2,...,m)
% f.et = et (et the parametrization of the boundary of the polygonal 
%            domain G) 
% f.etp = etp (etp is the first derivative of the parametrization et) 
% f.zet = zet (zet the parametrization of the boundary of the circular 
%              domain D; zet=f(et)) 
% f.zetp = zetp (zetp is the first derivative of the parametrization zet) 
% f.cent = cent (cent is a vector of length m where cent(k) is the center
%                of the circle C_k=f(\Gamma_k), k=1,2,...,m. 
%                If G is bounded, then cent(m)=0.)
% f.rad = rad (rad is a vector of length m where rad(k) is the radius
%                of the circle C_k=f(\Gamma_k), k=1,2,...,m. 
%                If G is bounded, then rad(m)=1.)
% f.imgver = imgver (imgver is a Cell Array where imgver{k}, k=1,2,...,m,
%                    contains the image of the vertices ver{k} under the 
%                    conformal map, i.e., imgver=f(ver)). 
% f.inf (only for unbounded G where f.inf=f'(inf)
% 
addpath pcm fmm % add the path of the required functions 
%
m            = length(ver); % number of polygons
n            = 2^9; % number of node points in each side of the polygons
                    % the accuracy of the results can be improved by 
                    % increasing the value of n. 
                    % However, choosing very large value of n could cause 
                    % a problem with the convergence (with the FMM method).
iprec        = 4; % the accuracy of the FMM method is 0.5e-12
gmresrestart = []; % the GMRES method is used without restart
gmrestol     = 0.5e-12; % tolerances of the GMRES method
gmresmaxit   = 100; % maximum number of iterations for the GMRES method
koebemaxit   = 100; % maximum number of iterations for the Koebe method
koebetol     = 1e-12; % tolerances of the Koebe method
%
if nargin == 2
    preimg1 = [];
end
% 
f = mainmap(ver,alpha,preimg1,m,n,iprec,gmresrestart,gmrestol,gmresmaxit,...
              koebetol,koebemaxit); % Compute the object f
%
end  