function   plotplri(f,n_of_cr,n_of_rd)
% plotplri.m 
% Nasser, September 6, 2019
%
% 
% This MATLAB function plots circular domain D with polar orthogonal grid
% and the polygonal domain G with the images of these polar grid under the 
% invers map, where
% n_of_cr: number of circles
% n_of_rd: number of rays
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
%%
if abs(alpha)<inf
    tang   =  linspace(0,2*pi,2000);
    trad   =  [linspace(0,0.9,1000),linspace(0.9,0.99999,1000)];
    for k=1:n_of_cr
        cirp{k} = (k/(n_of_cr+1)).*exp(i.*tang);
        for j=1:m-1
            cirp{k}(abs(cirp{k}-cent(j))./rad(j)<1+1e-6)=NaN+i*NaN;
        end
        cirpvc{k} = cirp{k}(abs(cirp{k})>=0);
        icirpvc{k} =  evalu(f,cirpvc{k},'v');
        icirp{k}  =NaN(size(cirp{k}))+i*NaN(size(cirp{k}));
        icirp{k}(abs(cirp{k})>=0)  =icirpvc{k};
    end
    for k=1:n_of_rd
        radp{k} =  trad.*exp(k*2*pi*i/n_of_rd);
        for j=1:m-1
            radp{k}(abs(radp{k}-cent(j))./rad(j)<1+1e-6)=NaN+i*NaN;
        end
        radpvc{k} = radp{k}(abs(radp{k})>=0);
        iradpvc{k} =  evalu(f,radpvc{k},'v');
        iradp{k}  =NaN(size(radp{k}))+i*NaN(size(radp{k}));
        iradp{k}(abs(radp{k})>=0)  =iradpvc{k};
    end
end
if abs(alpha)==inf
    mxr    =  max(real(zet));
    mir    =  min(real(zet));
    mxi    =  max(imag(zet));
    mii    =  min(imag(zet));
    neworg = (mxr+mir)/2+i*(mxi+mii)/2;
    mxrad  =  1.15*max(abs(zet-neworg)); 
    tang   =  linspace(0,2*pi,2000);
    trad   =  linspace(0,0.9999,2000);
    for k=1:n_of_cr
        cirp{k} = neworg+(mxrad*k/(n_of_cr)).*exp(i.*tang);
        for j=1:m
            cirp{k}(abs(cirp{k}-cent(j))./rad(j)<1+1e-6)=NaN+i*NaN;
        end
        cirpvc{k} = cirp{k}(abs(cirp{k})>=0);
        icirpvc{k} =  evalu(f,cirpvc{k},'v');
        icirp{k}  =NaN(size(cirp{k}))+i*NaN(size(cirp{k}));
        icirp{k}(abs(cirp{k})>=0)  =icirpvc{k};
    end    
    for k=1:n_of_rd
        radp{k} = neworg+mxrad*trad.*exp(k*2*pi*i/n_of_rd);
        for j=1:m
            radp{k}(abs(radp{k}-cent(j))./rad(j)<1+1e-6)=NaN+i*NaN;
        end
        radpvc{k} = radp{k}(abs(radp{k})>=0);
        iradpvc{k} =  evalu(f,radpvc{k},'v');
        iradp{k}  =NaN(size(radp{k}))+i*NaN(size(radp{k}));
        iradp{k}(abs(radp{k})>=0)  =iradpvc{k};
    end    
end
figure;
hold on
box on
axis equal
for k=1:length(ver)
    crv=et(1+sum(nv(1:k-1)):sum(nv(1:k)),1);
    plot(real(crv),imag(crv),'-k','LineWidth',1);
%     plot(real(ver{k}),imag(ver{k}),'sr','MarkerFaceColor','r');
end
for k=1:n_of_cr
    plot(real(icirp{k}),imag(icirp{k}),'-r');
end
for k=1:n_of_rd
    plot(real(iradp{k}),imag(iradp{k}),'-b');
end
% if abs(alpha)<inf
%     plot(real(alpha),imag(alpha),'or','MarkerFaceColor','r');
% end
if m==1
    ax_xp = max(real(et));
    ax_xm = min(real(et));
    ax_yp = max(imag(et));
    ax_ym = min(imag(et));
    if (ax_xp-ax_xm)>(ax_yp-ax_ym)
        LLx = 0.1*(ax_xp-ax_xm);
        LLy = ((ax_xp-ax_xm)-(ax_yp-ax_ym))/2+LLx;
        axis([ax_xm-LLx  ax_xp+LLx ax_ym-LLy  ax_yp+LLy])
    else
        LLy = 0.1*(ax_yp-ax_ym);
        LLx = ((ax_yp-ax_ym)-(ax_xp-ax_xm))/2+LLy;
        axis([ax_xm-LLx  ax_xp+LLx ax_ym-LLy  ax_yp+LLy])
    end
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
    plot(real(crv),imag(crv),'-k','LineWidth',1);
%     plot(real(imgver{k}),imag(imgver{k}),'sr','MarkerFaceColor','r');
end
for k=1:n_of_cr
    plot(real(cirp{k}),imag(cirp{k}),'-r');
end
for k=1:n_of_rd
    plot(real(radp{k}),imag(radp{k}),'-b');
end
% plot(real(cent(m)),imag(cent(m)),'or','MarkerfaceColor','r');
% plot(real(cent),imag(cent),'or','MarkerfaceColor','r');
% ax = gca;
% ax.XAxisLocation = 'origin';
% ax.YAxisLocation = 'origin';
%%
end