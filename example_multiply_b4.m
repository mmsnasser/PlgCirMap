% example_multiply_b4.m
% Nasser, September 5, 2019
% In this code, we consider:
% bounded multiply connected domain G of connectivity 4
% 
% 
clc
clear all
%
% The vertices of the polygon
% The inner polygons (the vertices must be clockwise oriented)
ver{1}=[ 1.5+1.0i ;  1.5+0.0i ;  0.5+0.0i ; 0.5+1.0i];
ver{2}=[-0.5+1.0i ; -1.5+0.0i ;-1.5+1.0i];
ver{3}=[-0.0-0.5i ;  0.5-1.0i ; 0.5-1.5i  ;-0.5-1.5i ;-0.5-1.0i];
% The outer polygon (the vertices must be counterclockwise oriented)
ver{4}=[ 2.0+2.0i ; -2.0+2.0i ; -2.0-2.0i ; 2.0-2.0i];
% Choose alpha, an auxiliary point in the domain G
alpha = 0;
%%
tic
% f=plgcirmap(ver,alpha);% f is the conformal mapping from the domain G
                         % onto the circular domain D with the  
                         % normalization f(alpha)=0 and f'(alpha)>0
f=plgcirmap(ver,alpha,ver{end}(3));% f is the conformal mapping from  
                         % the domain G onto the circular domain D with the  
                         % normalization f(alpha)=0 and f(ver{end}(k))=1
toc
%%
plotmap(f); % to plot the domain G and the circular domain D 
plotmap(f,'v','plr',21,21); % to plot polar grids in the circular domain
                              % D and their images in the domain G under  
                              % the invers map
plotmap(f,'v','rec',21,21); % to plot rectangular grids in the circular 
                              % domain D and their images in the domain G
                              % under the invers map
plotmap(f,'d','rec',21,21); % to plot rectangular grids in the domain 
                              % G and their images in the circular domain  
                              % D under the conformal map
plotmap(f,'d','plr',21,21); % to plot polar grids in the domain 
                              % G and their images in the circular domain  
                              % D under the conformal map
%%

%%
% Checking the accuracy of the toolbox PlgCirMap:
% 
% we choose test points: ztest in the domain G 
ztest   =  linspace(-1.999,1.999,999)-0.25i;
% ztest   =  linspace(-1.999,1.999,999)+1.9i;
% ztest   =  linspace(-1.999,1.999,999)+1.99i;
% ztest   =  linspace(-1.999,1.999,999)+1.999i;
% compute their images under the circular map: wtest
% We compute the images of the test points ztest under the conformal map 
% from G onto D 
wtest   =  evalu(f,ztest,'d');
% We compute the values of f^-1(f(zztest))  
ztest2  =  evalu(f,wtest,'v');
% We compute the maximum norm of the different between the test points 
% ztest and ztest2
error_norm = norm(ztest-ztest2,inf)
Error = abs(ztest-ztest2);
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
plot(real(ztest),imag(ztest),'.b')
%
figure;
hold on
box on
axis equal
for k=1:m
    crv=zet(1+sum(nv(1:k-1)):sum(nv(1:k)),1);
    plot(real(crv),imag(crv),'-k','LineWidth',2);
end
plot(real(wtest),imag(wtest),'.b')
%
figure
semilogy(real(ztest),Error,'b','LineWidth',1);
grid on
%%
