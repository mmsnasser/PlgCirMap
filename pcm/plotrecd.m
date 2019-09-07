function   plotrecd(f,n_of_hl,n_of_vl)
% plotrecd.m 
% Nasser, September 6, 2019
%
% 
% This MATLAB function plots polygon domain G with rectangular orthogonal 
% grid and the circular domain D with the images of these rectangular grid 
% under the dirct map, where
% n_of_hl: number of horizontal lines
% n_of_vl: number of vertical lines
%
%
nv    =  f.nv;
alpha =  f.alpha;
et    =  f.et;
zet   =  f.zet;
ver   =  f.ver;
cent  =  f.cent;
rad   =  f.rad;
m     =  length(ver);
%
%
%
%
if abs(alpha)<inf
    mxr    =  max(real(et));
    mir    =  min(real(et));
    mxi    =  max(imag(et));
    mii    =  min(imag(et));
    Lreal  =  mxr-mir;
    Limag  =  mxi-mii;
    tx     =  linspace(mir,mxr,1000);
    ty     =  linspace(mii,mxi,1000);
    for k=1:n_of_hl
        xpt{k} = tx+i*(Limag*k/(n_of_hl+1)+mii);
        for j=1:m-1
            inm = []; onm = [];
            [inm onm] = inpolygon(real(xpt{k}),imag(xpt{k}),real(ver{j}),imag(ver{j}));
            xpt{k}(inm)=NaN+i*NaN;  xpt{k}(onm) =NaN+i*NaN;
        end
        inm = []; onm = [];
        [inm onm] = inpolygon(real(xpt{k}),imag(xpt{k}),real(ver{m}),imag(ver{m}));
        xpt{k}(~inm)=NaN+i*NaN;  xpt{k}(onm) =NaN+i*NaN;
        xptvc{k}  = xpt{k}(abs(xpt{k})>=0);
        ixptvc{k} = evalu(f,xptvc{k},'d');
        ixpt{k}   = NaN(size(xpt{k}))+i*NaN(size(xpt{k}));
        ixpt{k}(abs(xpt{k})>=0)  =ixptvc{k};
    end
    for k=1:n_of_vl
        ypt{k} = Lreal*k/(n_of_vl+1)+mir+i*ty;
        for j=1:m-1
            inm = []; onm = [];
            [inm onm] = inpolygon(real(ypt{k}),imag(ypt{k}),real(ver{j}),imag(ver{j}));
            ypt{k}(inm)=NaN+i*NaN;  ypt{k}(onm) =NaN+i*NaN;
        end
        inm = []; onm = [];
        [inm onm] = inpolygon(real(ypt{k}),imag(ypt{k}),real(ver{m}),imag(ver{m}));
        ypt{k}(~inm)=NaN+i*NaN;  ypt{k}(onm) =NaN+i*NaN;
        yptvc{k}  = ypt{k}(abs(ypt{k})>=0);
        iyptvc{k} = evalu(f,yptvc{k},'d');
        iypt{k}  = NaN(size(ypt{k}))+i*NaN(size(ypt{k}));
        iypt{k}(abs(ypt{k})>=0)  =iyptvc{k};        
    end
end
if abs(alpha)==inf
    mxr    =  max(real(et));
    mir    =  min(real(et));
    mxi    =  max(imag(et));
    mii    =  min(imag(et));
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
    for k=1:n_of_hl
        xpt{k} = tx+i*(Limag*(k-1)/(n_of_hl-1)+mii);
        for j=1:m
            inm = []; onm = [];
            [inm onm] = inpolygon(real(xpt{k}),imag(xpt{k}),real(ver{j}),imag(ver{j}));
            xpt{k}(inm)=NaN+i*NaN;  xpt{k}(onm) =NaN+i*NaN;
        end
        xptvc{k}  = xpt{k}(abs(xpt{k})>=0);
        ixptvc{k} = evalu(f,xptvc{k},'d');
        ixpt{k}   = NaN(size(xpt{k}))+i*NaN(size(xpt{k}));
        ixpt{k}(abs(xpt{k})>=0)  =ixptvc{k};
    end
    for k=1:n_of_vl
        ypt{k} = Lreal*(k-1)/(n_of_vl-1)+mir+i*ty;
        for j=1:m
            inm = []; onm = [];
            [inm onm] = inpolygon(real(ypt{k}),imag(ypt{k}),real(ver{j}),imag(ver{j}));
            ypt{k}(inm)=NaN+i*NaN;  ypt{k}(onm) =NaN+i*NaN;
        end
        yptvc{k}  = ypt{k}(abs(ypt{k})>=0);
        iyptvc{k} = evalu(f,yptvc{k},'d');
        iypt{k}  = NaN(size(ypt{k}))+i*NaN(size(ypt{k}));
        iypt{k}(abs(ypt{k})>=0)  =iyptvc{k};        
    end
end
% 
% 
figure;
hold on
box on
axis equal
for k=1:length(ver)
    crv=et(1+sum(nv(1:k-1)):sum(nv(1:k)),1);
    plot(real(crv),imag(crv),'-k','LineWidth',1);
%     plot(real(ver{k}),imag(ver{k}),'sr','MarkerFaceColor','r');
end
for k=1:n_of_hl
    plot(real(xpt{k}),imag(xpt{k}),'-r');
end
for k=1:n_of_vl
    plot(real(ypt{k}),imag(ypt{k}),'-b');
end
% if abs(alpha)<inf
%     plot(real(alpha),imag(alpha),'or','MarkerFaceColor','r');
% end
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
% plot(real(cent(m)),imag(cent(m)),'or','MarkerfaceColor','r');
for k=1:n_of_hl
    plot(real(ixpt{k}),imag(ixpt{k}),'-r');
end
for k=1:n_of_vl
    plot(real(iypt{k}),imag(iypt{k}),'-b');
end
% ax = gca;
% ax.XAxisLocation = 'origin';
% ax.YAxisLocation = 'origin';
%%
end