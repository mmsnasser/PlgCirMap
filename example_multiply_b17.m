% example_multiply_b17.m
% Nasser, September 6, 2019
% In this code, we consider:
% bounded multiply connected domain G of connectivity 17
% 
% 
clc
clear all
%
% The vertices of the polygon
% The inner polygons (the vertices must be clockwise oriented)
ver{1}  = [31+10i ; 31+5i  ; 28+5i  ; 28+10i ];
ver{2}  = [25+10i ; 25+5i  ; 22+5i  ; 22+10i ];
ver{3}  = [19+10i ; 19+1i  ; 13+1i  ; 13+10i ];
ver{4}  = [10+10i ; 10+5i  ;  7+5i  ;  7+10i ];
ver{5}  = [ 4+10i ;  4+5i  ;  1+5i  ;  1+10i ];
ver{6}  = [31+19i ; 31+14i ; 28+14i ; 28+19i ];
ver{7}  = [25+19i ; 25+14i ; 22+14i ; 22+19i ];
ver{8}  = [19+14i ; 19+12i ; 17+12i ; 17+14i ];
ver{9}  = [15+14i ; 15+12i ; 13+12i ; 13+14i ];
ver{10} = [19+18i ; 19+16i ; 17+16i ; 17+18i ];
ver{11} = [15+18i ; 15+16i ; 13+16i ; 13+18i ];
ver{12} = [19+22i ; 19+20i ; 17+20i ; 17+22i ];
ver{13} = [15+22i ; 15+20i ; 13+20i ; 13+22i ];
ver{14} = [10+19i ; 10+14i ;  7+14i ;  7+19i ];
ver{15} = [ 4+19i ;  4+14i ;  1+14i ;  1+19i ];
ver{16} = [16+29i ; 23+24i ;  9+24i ];
% The outer polygon (the vertices must be counterclockwise oriented)
ver{17} = [16+32i ;  0+22i ;  0+0i  ; 32+0i  ; 32+22i ];
% 
% Choose alpha, an auxiliary point in the domain G
alpha = 16+15i;
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
plot(real(alpha),imag(alpha),'dr','MarkerfaceColor','r')
%
%
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
ztest   =  linspace(0.01,31.99,1000)+19.5i;
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

