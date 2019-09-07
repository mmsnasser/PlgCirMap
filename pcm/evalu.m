function ww = evalu (f,zz,b)
% evalu.m 
% Nasser, September 6, 2019
% 
% This MATLAB function compute the values ww of the analytic function f at
% the points zz in the polygon domain G.
% The function can compute also the values of the inverse function of f 
% at the points zz in the circular domain D
% 
% 
% Input:
%    1) f, an object with:
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
    z  = zz(:).';
    if abs(alpha)<inf
        w4  = fcaun(et,etp,zet,z,nv);
    elseif abs(alpha)==inf
        if isempty(f.inf)
            c2 = 1;
        else
            c2  = f.inf;
        end
        cnto(1) = getInteriorPoint(et(1:nv(1)), etp(1:nv(1)));
        w4  = (z-cnto(1)).*fcaun(et,etp,zet./(et-cnto(1)),z,nv,c2);
    end
    if isrow(zz)==1
        ww=w4(:).';
    else
        ww=w4(:);
    end
    %
elseif b=='v'
    %
    w4  = zz(:).';
    if abs(alpha)<inf
        z  = fcaun(zet,zetp,et,w4,nv);
    elseif abs(alpha)==inf
        if isempty(f.inf)
            c2v = 1;
        else
            c2v  = 1./f.inf;
        end
        z  = (w4-cent(1)).*fcaun(zet,zetp,et./(zet-cent(1)),w4,nv,c2v);
    end
    if isrow(zz)==1
        ww=z(:).';
    else
        ww=z(:);
    end
    %
elseif b=='db'
    %
    z   = zz(:).';
    ind = 0;
    zb  = [];
    for k4=1:length(z)
        [minj,j]=min(abs(z(k4)-et));
        if minj>=1e-6
            ind = ind+1;
            zb(ind)=z(k4);
        end
    end
    if length(zb)>0
        if abs(alpha)<inf
            wb  = fcaun(et,etp,zet,zb,nv);
        elseif abs(alpha)==inf
            if isempty(f.inf)
                c2 = 1;
            else
                c2  = f.inf;
            end
            cnto(1) = getInteriorPoint(et(1:nv(1)), etp(1:nv(1)));
            wb  = (zb-cnto(1)).*fcaun(et,etp,zet./(et-cnto(1)),zb,nv,c2);
        end
    end
    ind = 0;
    for k4=1:length(z)
        [minj,j]=min(abs(z(k4)-et));
        if minj<1e-6
            w4(k4,1)=zet(j);
        else
            ind = ind+1;
            w4(k4,1)=wb(ind);
        end
    end
    if isrow(zz)==1
        ww=w4(:).';
    else
        ww=w4(:);
    end
    %
elseif b=='vb'
    %
    w4  = zz(:).';
    ind = 0;
    wb  = [];
    for k4=1:length(w4)
        [minj,j]=min(abs(w4(k4)-zet));
        if minj>=1e-6
            ind = ind+1;
            wb(ind)=w4(k4);
        end
    end
    if length(wb)>0
        if abs(alpha)<inf
            zb  = fcaun(zet,zetp,et,wb,nv);
        elseif abs(alpha)==inf
            if isempty(f.inf)
                c2v = 1;
            else
                c2v  = 1./f.inf;
            end
            zb  = (wb-cent(1)).*fcaun(zet,zetp,et./(zet-cent(1)),wb,nv,c2v);
        end
    end
    ind = 0;
    for k4=1:length(w4)
        [minj,j]=min(abs(w4(k4)-zet));
        if minj<1e-6
            z(k4,1)=et(j);
        else
            ind = ind+1;
            z(k4,1)=zb(ind);
        end
    end
    if isrow(zz)==1
        ww=z(:).';
    else
        ww=z(:);
    end
    %
end
%
end
%%