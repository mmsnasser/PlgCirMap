function   f = mainmap(ver,alpha,preimg1,m,n,iprec,gmresrestart,gmrestol,gmresmaxit,koebetol,koebemaxit)
% mainmap.m 
% Nasser, September 6, 2019
%
% 
% This MATLAB function compute the conformal mapping w=f(z) from a polygon
% domain G onto a circular domin D and its invers z=f^-1(w). The domain G
% is multiply connected of connectivity m with the boundary
% \Gamma=\Gamma_1 U ... U \Gamma_m. If G is bounded, then the external
% boundary component is \Gamma_m.
% When m=1, the domain G is simply connected domain.
% 
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
% m  (number of polygons)
% n  (number of node points in each side of the polygons)
% iprec (the accuracy of the FMM method is 0.5e-12)
% gmresrestart = []; (the GMRES method is used without restart)
% gmrestol     (tolerances of the GMRES method)
% gmresmaxit   (maximum number of iterations for the GMRES method)
% koebemaxit   (maximum number of iterations for the Koebe method)
% koebetol     (tolerances of the Koebe method)
% 
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
% 
%
%
%
%
f.ver        =  ver;
f.alpha      =  alpha;
%
mm  =   2;
for k=1:m
    nv(k)    =   length(ver{k})*n;
    if k==m & abs(alpha)<inf
        nv(k)  =  mm*nv(k);
    end
end
for k=1:m
    Jk        =  1+sum(nv(1:k-1)):sum(nv(1:k));
    [et(Jk,1),etp(Jk,1)]=polygonp(ver{k},nv(k)/length(ver{k}));
end
%
f.nv    =  nv;
f.et    =  et;
f.etp   =  etp;
% 
if abs(alpha)<inf
    [zet,zetp,cent,rad] = cirmapb(et,etp,alpha,nv,iprec,gmresrestart,gmrestol,gmresmaxit,koebetol,koebemaxit);
elseif abs(alpha)==inf
    [zet,zetp,cent,rad] = cirmapu(et,etp,alpha,nv,iprec,gmresrestart,gmrestol,gmresmaxit,koebetol,koebemaxit);
    f.inf = [];
end
%
if (~isempty(preimg1) & abs(alpha)<inf)
    kk       =  length(ver{m});    
    bdzet    =  zet(sum(nv(1:m-1))+(kk-1)*mm*n+1);
    C2       =  exp(-i.*angle(bdzet));
    zet     =  C2.*zet;
    zetp    =  C2.*zetp;
    cent    =  C2.*cent;
end
% 
if (~isempty(preimg1) & abs(alpha)==inf)
    kk       =  length(ver{m});    
    bdzet    =  zet(sum(nv(1:m-1))+(kk-1)*n+1);
    C1       =  cent(m);
    C2       =  exp(-i.*angle(bdzet-cent(m)))./rad(m);
    zet      =  C2.*(zet -C1);
    zetp     =  C2.*(zetp);
    cent     =  C2.*(cent-C1);
    rad      =  rad./rad(m);
    f.inf  =  C2;
end
%
f.zet      =  zet;
f.zetp     =  zetp;
f.cent     =  cent;
f.rad      =  rad;
%
if abs(alpha)<inf
    for k=1:m-1
        for j=1:length(ver{k})
            imgver{k}(j,1) =  zet(sum(nv(1:k-1))+(j-1)*n+1);
        end
    end
    for j=1:length(ver{m}) %k=m
        imgver{m}(j,1) =  zet(sum(nv(1:m-1))+(j-1)*mm*n+1);
    end
elseif abs(alpha)==inf
    for k=1:m
        for j=1:length(ver{k})
            imgver{k}(j,1) =  zet(sum(nv(1:k-1))+(j-1)*n+1);
        end
    end
end
%
f.imgver   =  imgver;
%%
end  