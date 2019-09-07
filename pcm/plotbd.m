function   plotbd(f)
% plotbd.m 
% Nasser, September 6, 2019
%
% 
% This MATLAB function plots polygon domain G and the circular domain D;
% and the vertices of the polygons with their images on the circles
%
%
nv     =  f.nv;
et     =  f.et;
zet    =  f.zet;
cent   =  f.cent;
alpha  =  f.alpha;
imgver =  f.imgver;
rad    =  f.rad;
ver    =  f.ver;
m      =  length(ver);
%
%
figure;
colormap prism;
cmap=colormap;
hold on
box on
axis equal
for k=1:length(ver)
    crv=et(1+sum(nv(1:k-1)):sum(nv(1:k)),1);
    plot(real(crv),imag(crv),'-k','LineWidth',2);
    for j=1:length(ver{k})
        jj = (1:length(ver{k}))';
        plot(real(ver{k}(j)),imag(ver{k}(j)),'sk','MarkerFaceColor',cmap(1*j,:));
    end
end
if abs(alpha)<inf
    plot(real(alpha),imag(alpha),'pr','MarkerFaceColor','r');
    text(real(alpha)+0.05,imag(alpha)+0.1,'{$\alpha$}','FontSize',14,'Interpreter','latex');
end
% ax = gca;
% ax.XAxisLocation = 'origin';
% ax.YAxisLocation = 'origin';
%
%
figure;
hold on
box on
axis equal
for k=1:m
    crv=zet(1+sum(nv(1:k-1)):sum(nv(1:k)),1);
    plot(real(crv),imag(crv),'-k','LineWidth',2);
%     plot(real(cent(k)),imag(cent(k)),'o','LineWidth',2);
    for j=1:length(ver{k})
        jj = (1:length(ver{k}))';
        plot(real(imgver{k}(j)),imag(imgver{k}(j)),'sk','MarkerFaceColor',cmap(1*j,:));
    end
end
% ax = gca;
% ax.XAxisLocation = 'origin';
% ax.YAxisLocation = 'origin';
%
%
end