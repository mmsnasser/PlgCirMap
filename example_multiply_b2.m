% example_multiply_b2.m
% Nasser, September 5, 2019
% In this code, we consider:
% bounded multiply connected domain G of connectivity 2
% 
% 
clc
clear all
%
% The vertices of the polygon
% The inner polygon (the vertices must be clockwise oriented)
ver{1}=[ 1.5+1.0i ;  1.5      ;  0.5      ; 0.5+1.0i];
% The outer polygon (the vertices must be counterclockwise oriented)
ver{2}=[ 2.0+2.0i ; -3.0+2.0i ; -2.0-2.0i ; 2.0-2.0i];
% Choose alpha, an auxiliary point in the domain G
alpha = 0;
%%
tic
f=plgcirmap(ver,alpha);% f is the conformal mapping from the domain G
                         % onto the circular domain D with the  
                         % normalization f(alpha)=0 and f'(alpha)>0
% f=plgcirmap(ver,alpha,ver{end}(end));% f is the conformal mapping from  
                         % the domain G onto the circular domain D with the  
                         % normalization f(alpha)=0 and f(ver{end}(end))=1
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
ztest  =  0.0+i.*linspace(-1.99,1.99,1000);
% We compute the images of the test points ztest under the conformal map 
% from G onto D 
wtest  =  evalu(f,ztest,'d');
% We compute the values of f^-1(f(zztest))  
ztest2 =  evalu(f,wtest,'v');
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
semilogy(imag(ztest),Error,'b','LineWidth',1);
grid on
%%
% Checking the accuracy of the toolbox PlgCirMap:
% 
% we choose test points: wtesti in the circular domain D 
wtesti  =  exp(-i*0.25*pi).*linspace(-0.99,0.99,1000)+0*i;
% We compute the images of the test points wtesti under the inverse  
% conformal map from D onto G
ztesti  =  evalu(f,wtesti,'v');
% We compute the values of f(f^-1(ztesti))
wtesti2 =  evalu(f,ztesti,'d');
% We compute the maximum norm of the different between the test points 
% wtesti and wtesti2
error_norm = norm(wtesti-wtesti2,inf)
Errori = abs(wtesti-wtesti2);
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
plot(real(ztesti),imag(ztesti),'.b')
%
figure;
hold on
box on
axis equal
for k=1:m
    crv=zet(1+sum(nv(1:k-1)):sum(nv(1:k)),1);
    plot(real(crv),imag(crv),'-k','LineWidth',2);
end
plot(real(wtesti),imag(wtesti),'.b')
%
figure
semilogy(real(wtesti),Errori,'b','LineWidth',1);
grid on
%%