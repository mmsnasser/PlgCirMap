function   [et,etp]=polygonp(ver,ns)
% polygonp.m
% Nasser, June 10, 2019
% This function compute the discretization of the parametrization of the
% polygon with the vertices "ver' where "ns" is the graded points on each
% side of the polygon. The number of graded points for the whole polygon
% is "n=number os sides*ns". The graded mesh points are computed by the 
% function "deltw.m".
% 
% 
n      =   length(ver)*ns;
t      =  (0:2*pi/n:2*pi-2*pi/n).';
[s,sp] =   deltw(t,length(ver),3);
for j=1:length(ver)
    sv{j}  =  s((j-1)*n/length(ver)+1:j*n/length(ver));
end
verts     =   ver; verts(length(ver)+1)=verts(1);
for j=1:length(ver)
    etv{j}   = verts(j)+(length(ver)/(2*pi))*(verts(j+1)-verts(j)).*(sv{j}-sv{j}(1));
    etvp{j}  =          (length(ver)/(2*pi))*(verts(j+1)-verts(j)).*(ones(size(sv{j})));
end
eto = []; etopo = [];
for j=1:length(ver)
    eto((j-1)*n/length(ver)+1:j*n/length(ver),1)     =  etv{j};
    etopo((j-1)*n/length(ver)+1:j*n/length(ver),1)   =  etvp{j};
end
et  =  eto; etp =  etopo.*sp;
end
%