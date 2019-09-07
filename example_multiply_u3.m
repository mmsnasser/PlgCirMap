% example_multiply_u3.m
% Nasser, September 6, 2019
% In this code, we consider:
% unbounded multiply connected domain G of connectivity 3
% 
% 
clc
clear all
%
% The vertices of the polygons (the vertices must be clockwise oriented)
ver{1}=[ 1.5+1.0i ;  1.5+0.0i ;  0.5+0.0i ; 0.5+1.0i];
ver{2}=[-0.5+1.0i ; -1.5+0.0i ;-1.5+1.0i];
ver{3}=[-0.0-0.5i ;  0.5-1.0i ; 0.5-1.5i  ;-0.5-1.5i ;-0.5-1.0i];
% The domain G is unbounded. So, we choose alpha=inf.
alpha = inf;
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
ttest   =  linspace(0,2*pi,1000);
ztest   =  2.*exp(i.*ttest);
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
semilogy(ttest,Error,'b','LineWidth',1);
grid on
%%