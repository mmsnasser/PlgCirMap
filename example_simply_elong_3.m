% example_simply_enlog_2.m
% Nasser, September 3, 2019
% In this code, we consider:
% bounded simply connected domain G inside elongated polygon with 8 vertices
clc
clear all
%
% The vertices of the polygon
ver{1} = [ 6+1i; -4+1i; -4+9i; 6+9i; 6+11i; -6+11i ; -6-1i ;  6-1i];
% Choose alpha, an auxiliary point in the domain G
alpha =  -5.0+5.0i;
%%
tic
f=plgcirmap(ver,alpha);% f is the conformal mapping from the domain G
                         % onto the circular domain D with the  
                         % normalization f(alpha)=0 and f'(alpha)>0
% f=plgcirmap(ver,alpha,ver{1}(end));% f is the conformal mapping from the 
                         % domain G onto the circular domain D with the  
                         % normalization f(alpha)=0 and f(ver{1}(end))=1
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
% Checking the accuracy of the toolbox PlgCirMap:
% 
% we choose test points: wtest in the circular domain D 
ttest = linspace(0,2*pi,1000);
wtest = 0.9.*exp(i.*ttest);
% We compute the values of the test points wtest under the inverse map from
% D onto G using the toolbox PlgCirMap
ztestie = evalu(f,wtest,'v');
% We compute the maximum norm of the different between the test points 
% wtest and the computed values of f(f^-(wtest)) for PlgCirMap toolbox
error_ie = norm(wtest-evalu(f,ztestie,'d'),inf)
%%
% Checking the accuracy of the toolbox PlgCirMap:
% 
% we choose test points: zztest in the circular domain G 
zztest  = alpha+0.4.*exp(i.*ttest);
% We compute the values of the test points zztest under the conformal map 
% from G onto D using the toolbox PlgCirMap
wwtestie = evalu(f,zztest,'d');
% We compute the maximum norm of the different between the test points 
% zztest and the computed values of f^-1(f(zztest)) for PlgCirMap toolbox 
error_ie = norm(zztest-evalu(f,wwtestie,'v'),inf)
%%