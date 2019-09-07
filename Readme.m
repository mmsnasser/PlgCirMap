function Readme()
%README
%
% Mohamed Nasser, 2019
% Please cite these collections of MATLAB files as:

% M.M.S. Nasser, PlgCirMap: A MATLAB toolbox for computing the conformal 
% maps from polygonal multiply connected domains onto circular domains. 
% https://github.com/mmsnasser/plgcirmap.
%
% 
% The main file in these collections of MATLAB files is the file
% plgcirmap.m
% 
%    f = plgcirmap(ver,alpha,preimg1)
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
% In these collection of MATLAB files, we use the MATLAB function 
% fbie.m
% which is available in:
% 
% M.M.S. Nasser, FBIEGNK:  A MATLAB toolbox for fast solution of 
% boundary integral equations with the generalized Neumann kernel, 
% Version 1.1, 2016. https://github.com/mmsnasser/FBIEGNK.
% 
% and in:
% M.M.S. Nasser, Fast solution of boundary integral equations with the 
% generalized Neumann kernel, Electronic Transactions on Numerical 
% Analysis,  44 (2015) 189--229.
% 
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% The MATLAB function fbie and other functions in this toolbox require
% the following files:
% 
% zfmm2dpart.m
% fmm2d_r2012a.mexw32
% fmm2d_r2012a.mexw64
% pthreadGC2-w32.dll
% pthreadGC2-w64.dll
% 
% from the MATLAB Toolbox:
% L. GREENGARD AND Z. GIMBUTAS , FMMLIB2D: A MATLAB toolbox for
% fast multipole method in two dimensions, Version 1.2, 2012.
%
% You can download the whole toolbox from:
% http://www.cims.nyu.edu/cmcl/fmm2dlib/fmm2dlib.html
% or from
% https://github.com/zgimbutas/fmmlib2d
%
% Please see: 
% https://github.com/zgimbutas/fmmlib2d/blob/master/COPYING
% for more details.
% 
% PLEASE cite the FMMLIB2D toolbox whenever you use the PlgCirMap toolbox
%
% Acknowledgments:
% I would like to thank Prof. Leslie Greengard and Prof. Zydrunas Gimbutas 
% for making the MATLAB toolbox FMMLIB2D publicly available.
%
%
% This program is free software; you can redistribute it and/or modify 
% it under the terms of the GNU General Public License as published by 
% the Free Software Foundation; either version 2 of the License, or 
% (at your option) any later version.  This program is distributed in 
% the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
% even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
% PARTICULAR PURPOSE.  See the GNU General Public License for more 
% details. You should have received a copy of the GNU General Public 
% License along with this program; 
% if not, see <http://www.gnu.org/licenses/>.
%
end