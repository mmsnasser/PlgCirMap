function   plotreci(f,n_of_h1,n_of_v1)
% plotreci.m 
% Nasser, September 6, 2019
%
% 
% This MATLAB function plots circular domain D with rectangular orthogonal 
% grid and the polygonal domain G with the images of these rectangular 
% grid under the invers map, where
% n_of_hl: number of horizontal lines
% n_of_vl: number of vertical lines
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
    tx     =  linspace(-1,1,1000);
    ty     =  linspace(-1,1,1000);    
    for k=1:n_of_h1
        hlin{k} = (2*k/(n_of_h1+1)-1)*i+tx;
        hlin{k}(abs(hlin{k})>1-1e-6)=NaN+i*NaN;
        for j=1:m-1
            hlin{k}(abs(hlin{k}-cent(j))./rad(j)<1+1e-6)=NaN+i*NaN;
        end
        hlinvc{k} = hlin{k}(abs(hlin{k})>=0);
        ihlinvc{k} =  evalu(f,hlinvc{k},'v');
        ihlin{k}  =NaN(size(hlin{k}))+i*NaN(size(hlin{k}));
        ihlin{k}(abs(hlin{k})>=0)  =ihlinvc{k};
    end
    for k=1:n_of_v1
        vlin{k} = (2*k/(n_of_v1+1)-1)+i*ty;
        vlin{k}(abs(vlin{k})>1-1e-6)=NaN+i*NaN;
        for j=1:m-1
            vlin{k}(abs(vlin{k}-cent(j))./rad(j)<1+1e-6)=NaN+i*NaN;
        end
        vlinvc{k} = vlin{k}(abs(vlin{k})>=0);
        ivlinvc{k} =  evalu(f,vlinvc{k},'v');
        ivlin{k}  =NaN(size(vlin{k}))+i*NaN(size(vlin{k}));
        ivlin{k}(abs(vlin{k})>=0)  =ivlinvc{k};
    end
end
if abs(alpha)==inf
    mxr    =  max(real(zet));
    mir    =  min(real(zet));
    mxi    =  max(imag(zet));
    mii    =  min(imag(zet));
    Lreal  =  mxr-mir;
    Limag  =  mxi-mii;
    if m==1
        Lreal=2*Lreal;
        Limag=2*Limag;
    end
    mxr    =  mxr+0.15*Lreal;
    mir    =  mir-0.15*Lreal;
    mxi    =  mxi+0.15*Limag;
    mii    =  mii-0.15*Limag;    
    Lreal  =  mxr-mir;
    Limag  =  mxi-mii;
    tx     =  linspace(mir,mxr,1000);
    ty     =  linspace(mii,mxi,1000);    
    for k=1:n_of_h1
        hlin{k} = tx+i*(Limag*(k-1)/(n_of_h1-1)+mii);
        for j=1:m
            hlin{k}(abs(hlin{k}-cent(j))./rad(j)<1+1e-6)=NaN+i*NaN;
        end
        hlinvc{k} = hlin{k}(abs(hlin{k})>=0);
        ihlinvc{k} =  evalu(f,hlinvc{k},'v');
        ihlin{k}  =NaN(size(hlin{k}))+i*NaN(size(hlin{k}));
        ihlin{k}(abs(hlin{k})>=0)  =ihlinvc{k};
    end    
    for k=1:n_of_v1
        vlin{k} = Lreal*(k-1)/(n_of_v1-1)+mir+i*ty;
        for j=1:m
            vlin{k}(abs(vlin{k}-cent(j))./rad(j)<1+1e-6)=NaN+i*NaN;
        end
        vlinvc{k} = vlin{k}(abs(vlin{k})>=0);
        ivlinvc{k} =  evalu(f,vlinvc{k},'v');
        ivlin{k}  =NaN(size(vlin{k}))+i*NaN(size(vlin{k}));
        ivlin{k}(abs(vlin{k})>=0)  =ivlinvc{k};
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
for k=1:n_of_h1
    plot(real(ihlin{k}),imag(ihlin{k}),'-r');
end
for k=1:n_of_v1
    plot(real(ivlin{k}),imag(ivlin{k}),'-b');
end
% if abs(alpha)<inf
%     plot(real(alpha),imag(alpha),'or','MarkerFaceColor','r');
% end
if m==1 & abs(alpha)<inf
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
for k=1:n_of_h1
    plot(real(hlin{k}),imag(hlin{k}),'-r');
end
for k=1:n_of_v1
    plot(real(vlin{k}),imag(vlin{k}),'-b');
end
% plot(real(cent(m)),imag(cent(m)),'or','MarkerfaceColor','r');
% plot(real(cent),imag(cent),'or','MarkerfaceColor','r');
% ax = gca;
% ax.XAxisLocation = 'origin';
% ax.YAxisLocation = 'origin';
%%
end