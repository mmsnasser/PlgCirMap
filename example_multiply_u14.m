% example_multiply_u14.m
% Nasser, September 6, 2019
% In this code, we consider:
% unbounded multiply connected domain G of connectivity 4
% 
% 
clc
clear all
%
% The vertices of the polygons (the vertices must be clockwise oriented)
ver{1}  = [ 4+4i ; 4+2i ; 2+2i ];
ver{2}  = [ 4+1i ; 4+0i ; 2+0i ; 2+1i ];
ver{3}  = [ 3-1i ; 3-2i ; 2-2i ; 2-1i ];
ver{4}  = [ 4-2i ; 4-4i ; 2-4i ];
ver{5}  = [ 2+3i ; 1+2i ; 0+3i ; 0+4i ; 1+4i ];
ver{6}  = [ 1+1i ; 0+0i ; 0+2i ];
ver{7}  = [ 1+0i ; 1-2i ; 0-2i ; 0-1i ];
ver{8}  = [ 2-3i ; 0-4i ; 0-3i ];
ver{9}  = [-1+4i ;-1+1i ;-2+1i ];
ver{10} = [-1+0i ;-1-1i ;-2-1i ;-2+0i ];
ver{11} = [-1-2i ;-1-4i ;-2-4i ;-2-3i ];
ver{12} = [-2+3i ;-2+2i ;-3+1i ;-4+2i ;-4+3i ;-3+4i ];
ver{13} = [-3+0i ;-2-2i ;-4-2i ;-4+0i ];
ver{14} = [-3-3i ;-3-4i ;-4-4i ;-4-3i ];
% 
% The domain G is unbounded. So, we choose alpha=inf.
alpha = inf;
%
% 
% Simple check of the domain to make sure that there is no overlap between
% the polygons and that the orientations of the polygons are clockwise
figure
hold on
axis equal
for k=1:length(ver)
    plgk = []; plgk = ver{k}; plgk(end+1)=plgk(1);
    plot(real(plgk),imag(plgk),'b')
    plot(real(ver{k}),imag(ver{k}),'-or')
end
%
%
%%
tic
f=plgcirmap(ver,alpha);% f is the conformal mapping from the domain G
                       % onto the circular domain D with the  
                       % normalization f(z)=z+O(1/z) near infinity
% f=plgcirmap(ver,alpha,ver{end}(end));% f is the conformal mapping from  
                       % the domain G onto the circular domain D with the  
                       % normalization: f(inf)=inf, cent(m)=0, rad(m)=1, 
                       % and f(ver{end}(end))=1
toc
%%
plotmap(f); % to plot the domain G and the circular domain D 
%%
plotmap(f,'v','plr'); % to plot polar grids in the circular domain
                              % D and their images in the domain G under  
                              % the invers map
%%                              
plotmap(f,'v','rec',21,21); % to plot rectangular grids in the circular 
                              % domain D and their images in the domain G
                              % under the invers map
%%                              
plotmap(f,'d','rec'); % to plot rectangular grids in the domain 
                              % G and their images in the circular domain  
                              % D under the conformal map
%%                              
plotmap(f,'d','plr',21,21); % to plot polar grids in the domain 
                              % G and their images in the circular domain  
                              % D under the conformal map
%%
% Checking the accuracy of the toolbox PlgCirMap:
% 
% we choose test points: ztest in the domain G 
tt      =  linspace(0,2*pi,1000);
zz      =  5.8.*exp(i.*tt);
% We compute the images of the test points ztest under the conformal map 
% from G onto D 
wzz     =  evalu(f,zz,'d');
% We  compute the values of f^-1(f(zztest))  
zzi     =  evalu(f,wzz,'v');
% We compute the maximum norm of the different between the test points 
% ztest and ztest2
error_norm = norm(zz-zzi,inf)
Error = abs(zz-zzi);
% 
nv     =  f.nv;
et     =  f.et;
zet    =  f.zet;
imgver =  f.imgver;
m      =  length(ver);
% 
figure;
hold on
box on
axis equal
for k=1:m
    crv=et(1+sum(nv(1:k-1)):sum(nv(1:k)),1);
    plot(real(crv),imag(crv),'-k','LineWidth',2);
end
plot(real(zz),imag(zz),'.b')
%
figure;
hold on
box on
axis equal
for k=1:m
    crv=zet(1+sum(nv(1:k-1)):sum(nv(1:k)),1);
    plot(real(crv),imag(crv),'-k','LineWidth',2);
end
plot(real(wzz),imag(wzz),'.b')
%
figure
semilogy(tt,Error,'b','LineWidth',1);
axis([0 2*pi 1e-9 1e-6])
grid on
%%
% error for the difference between the distnace from the boundary of the
% circles to their centers and the radius of the circles
for k=1:m
    err_centk(k)=norm(abs(zet(1+sum(nv(1:k-1)):sum(nv(1:k)))-f.cent(k))-f.rad(k),inf);
end
err_cent = norm(err_centk,inf)
%%