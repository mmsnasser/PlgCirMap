function w = evalud (f,z,b)
% evalud.m 
% Nasser, September 6, 2019
% 
% Let f be the conformal mapping from the polygonal domain G onto the
% circular domain D and g be its invers.
% For b='d', this function computes w=f'(z) where z is a vector of points
% in G 
% For b='v', this function computes w=g'(z) where z is a vector of points
% in D 
% 
% 
% Input:
%  1) f, an object with:
% 
% f.alpha = alpha (alpha is an auxiliary point in the domain G if G is
%                  bounded and alpha=inf if G is unbounded) 
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
% f.inf (only for unbounded G where f.inf=f'(inf)
% 
% 2) zz (zz is a vector of points. zz can be in the poygon domain G to 
%        compute their images under the dirct map; or zz can be in the 
%        circular domain D to compute their images under the inverse map.
%        
%        This function can also compute the values of the direct map and
%        the inverse map when the points zz are on the boundary. However,
%        the method used for such case need to be improved. So, this 
%        function will be improved in the next version to compute the 
%        values of the direct or the invers maps if the poits zz are on the 
%        boundary of the polygon domain G or the boundary of the circular 
%        domain D.
% 
% 3) b (b='d' for computing the direct map when zz in G;
%       b='v' for computing the inverse map when zz in D;
%       b='db' for computing the direct map when zz in the boundary of G;
%       b='vb' for computing the inverse map when zz in the boundary of D).

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
% 
% Output: ww (the computed values of the function f or its inverse at the
%             points zz).
% 
%
nv     =  f.nv;
alpha  =  f.alpha;
et     =  f.et;
etp    =  f.etp;
zet    =  f.zet;
zetp   =  f.zetp;
cent   =  f.cent;
%
%
if b=='d'
    %
    zz  = z(:).';
    if abs(alpha)<inf
        w4  = fcaunp(et,etp,zetp,zz,nv);
    elseif abs(alpha)==inf
        if isempty(f.inf)
            c2 = 1;
        else
            c2  = f.inf;
        end
        w4  = fcaunp(et,etp,zetp,zz,nv,c2);
    end
    if isrow(z)==1
        w=w4(:).';
    else
        w=w4(:);
    end
    %
elseif b=='v'
    %
    w4  = z(:).';
    if abs(alpha)<inf
        zz  = fcaunp(zet,zetp,etp,w4,nv);
    elseif abs(alpha)==inf
        if isempty(f.inf)
            c2v = 1;
        else
            c2v  = 1./f.inf;
        end
        zz  = fcaunp(zet,zetp,etp,w4,nv,c2v);
    end
    if isrow(z)==1
        w=zz(:).';
    else
        w=zz(:);
    end
    %
end
%
end
%%