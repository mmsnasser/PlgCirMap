% compare_sc_example_simply_enlog_2.m
% Nasser, September 3, 2019
% In this code, we consider:
% bounded simply connected domain G inside elongated polygon with 6 vertices
% 
% For comparison with Schwarz–Christoffel Toolbox, you need to download 
% the SC toolbox from http://www.math.udel.edu/~driscoll/SC/, then extract
% the toolbox into a folder sc, then add this folder to the folder contains
% this MATLAB file.
% 
clc
clear all
%
% The vertices of the polygon
ver{1} = [ 6+1i ; -4+1i ; -4+11i ; -6+11i ; -6-1i ;  6-1i];
% Choose alpha, an auxiliary point in the domain G
alpha =  -5.0+0.0i;
%% 
tic
% f=plgcirmap(ver,alpha);% f is the conformal mapping from the domain G
                         % onto the circular domain D with the  
                         % normalization f(alpha)=0 and f'(alpha)>0
f=plgcirmap(ver,alpha,ver{1}(end));% f is the conformal mapping from the 
                         % domain G onto the circular domain D with the  
                         % normalization f(alpha)=0 and f(ver{1}(end))=1
toc
%%
% plotmap(f); % to plot the domain G and the circular domain D 
% plotmap(f,'v','plr',21,21); % to plot polar grids in the circular domain
                              % D and their images in the domain G under  
                              % the invers map
% plotmap(f,'v','rec',21,21); % to plot rectangular grids in the circular 
                              % domain D and their images in the domain G
                              % under the invers map
% plotmap(f,'d','rec',21,21); % to plot rectangular grids in the domain 
                              % G and their images in the circular domain  
                              % D under the conformal map
% plotmap(f,'d','plr',21,21); % to plot polar grids in the domain 
                              % G and their images in the circular domain  
                              % D under the conformal map
%%
addpath sc % to use the Schwarz–Christoffel Toolbox 
options = scmapopt('Tolerance',1e-14); % the accuracet for SC toolbox
p=polygon(ver{1}); 
tic
fsc = diskmap(p,options);
toc
fsc = center(fsc,alpha);% fsc=f^-1 is the invers map from the circular 
                        % domain D onto the domain G compute by SC toolbox
                        % such that f(alpha)=0, f(ver{1}(end))=1
%%
% to plot polar grids in the circular domain D and their images in the 
% domain G under the invers map using both methods
plotmap(f,'v','plr',21,21);
figure;
plot(fsc,21,21)
%%
% Checking the accuracy of the toolbox PlgCirMap:
% 
% we choose test points: wtest in the circular domain D 
ttest = linspace(0,2*pi,1000);
wtest = 0.9.*exp(i.*ttest);
% We compute the values of the test points wtest under the inverse map from
% D onto G using the toolbox PlgCirMap
ztestie = evalu(f,wtest,'v');
% and using SC toolbox
ztestsc = fsc(wtest);
% We compute the maximum norm of the different between the computed values
error_inv_map = norm(ztestie-ztestsc,inf)
% We compute the maximum norm of the different between the test points 
% wtest and the computed values of f(f^-(wtest)) for SC toolbox
error_sc = norm(wtest-eval(inv(fsc),ztestsc),inf)
% We compute the maximum norm of the different between the test points 
% wtest and the computed values of f(f^-(wtest)) for PlgCirMap toolbox
error_ie = norm(wtest-evalu(f,ztestie,'d'),inf)
%%
% Checking the accuracy of the toolbox PlgCirMap:
% 
% we choose test points: zztest in the circular domain G 
zztest  = 0.4.*exp(i.*ttest);
% We compute the values of the test points zztest under the conformal map 
% from G onto D using the toolbox PlgCirMap
wwtestie = evalu(f,zztest,'d');
% We compute the values of the test points zztest under the conformal map 
% from G onto D using the SC toolbox 
wwtestsc = eval(inv(fsc),zztest);
% We compute the maximum norm of the different between the computed values
error_map = norm(wwtestie-wwtestsc,inf)
% We compute the maximum norm of the different between the test points 
% zztest and the computed values of f^-1(f(zztest)) for SC toolbox
error_sc = norm(zztest-fsc(wwtestsc),inf)
% We compute the maximum norm of the different between the test points 
% zztest and the computed values of f^-1(f(zztest)) for PlgCirMap toolbox 
error_ie = norm(zztest-evalu(f,wwtestie,'v'),inf)
%%
% Checking the accuracy of the toolbox PlgCirMap:
% 
% The prevertices computed using SC toolbox
prevertsc = get(fsc,'prevert');
% The prevertices computed using PlgCirMap toolbox
prevertie = f.imgver{1};
% We compute the maximum norm of the different between the computed values
error_prevert = norm(prevertsc-prevertie,inf)
%%