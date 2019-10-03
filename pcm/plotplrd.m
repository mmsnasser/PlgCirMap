function   plotplrd(f,n_of_cr,n_of_rd)
% plotplrd.m 
% Nasser, September 6, 2019
%
% 
% This MATLAB function plots polygon domain G with polar orthogonal grid
% and the circular domain D with the images of these polar grid under the 
% the dirct map, where
% n_of_cr: number of circles
% n_of_rd: number of rays
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
    neworg = (mxr+mir)/2+i*(mxi+mii)/2;
    mxrad  =  sqrt(((mxr-mir)/2)^2+((mxi-mii)/2)^2); 
    tang   =  linspace(0,2*pi,1000);
    trad   =  linspace(0,0.9999,1000);
    for k=1:n_of_cr
        cirp{k} = neworg+mxrad*(k/(n_of_cr)).*exp(i.*tang);
        for j=1:m-1
            inm = []; onm = [];
            [inm onm] = inpolygon(real(cirp{k}),imag(cirp{k}),real(ver{j}),imag(ver{j}));
            cirp{k}(inm)=NaN+i*NaN;  cirp{k}(onm) =NaN+i*NaN;
        end
        inm = []; onm = [];
        [inm onm] = inpolygon(real(cirp{k}),imag(cirp{k}),real(ver{m}),imag(ver{m}));
        cirp{k}(~inm)=NaN+i*NaN;  cirp{k}(onm) =NaN+i*NaN;
        cirpvc{k}  = cirp{k}(abs(cirp{k})>=0);
        icirpvc{k} = evalu(f,cirpvc{k},'d');
        icirp{k}   = NaN(size(cirp{k}))+i*NaN(size(cirp{k}));
        icirp{k}(abs(cirp{k})>=0)  =icirpvc{k};
    end
    for k=1:n_of_rd
        radp{k} =  neworg+mxrad*trad.*exp(k*2*pi*i/n_of_rd);
        for j=1:m-1
            inm = []; onm = [];
            [inm onm] = inpolygon(real(radp{k}),imag(radp{k}),real(ver{j}),imag(ver{j}));
            radp{k}(inm)=NaN+i*NaN;  radp{k}(onm) =NaN+i*NaN;
        end
        inm = []; onm = [];
        [inm onm] = inpolygon(real(radp{k}),imag(radp{k}),real(ver{m}),imag(ver{m}));
        radp{k}(~inm)=NaN+i*NaN;  radp{k}(onm) =NaN+i*NaN;
        radpvc{k}  = radp{k}(abs(radp{k})>=0);
        iradpvc{k} = evalu(f,radpvc{k},'d');
        iradp{k}  = NaN(size(radp{k}))+i*NaN(size(radp{k}));
        iradp{k}(abs(radp{k})>=0)  =iradpvc{k};        
    end
end
if abs(alpha)==inf
    mxr    =  max(real(et));
    mir    =  min(real(et));
    mxi    =  max(imag(et));
    mii    =  min(imag(et));
    neworg = (mxr+mir)/2+i*(mxi+mii)/2;
    mxrad  =  1.15*max(abs(et-neworg)); 
    if m==1
        mxrad=1.3*mxrad;
    end
    tang   =  linspace(0,2*pi,1000);
    trad   =  linspace(0,0.9999,1000);
    for k=1:n_of_cr
        cirp{k} = neworg+mxrad*(k/(n_of_cr)).*exp(i.*tang);
        for j=1:m
            inm = []; onm = [];
            [inm onm] = inpolygon(real(cirp{k}),imag(cirp{k}),real(ver{j}),imag(ver{j}));
            cirp{k}(inm)=NaN+i*NaN;  cirp{k}(onm) =NaN+i*NaN;
        end
        cirpvc{k}  = cirp{k}(abs(cirp{k})>=0);
        icirpvc{k} = evalu(f,cirpvc{k},'d');
        icirp{k}   = NaN(size(cirp{k}))+i*NaN(size(cirp{k}));
        icirp{k}(abs(cirp{k})>=0)  =icirpvc{k};
    end
    for k=1:n_of_rd
        radp{k} =  neworg+mxrad*trad.*exp(k*2*pi*i/n_of_rd);
        for j=1:m
            inm = []; onm = [];
            [inm onm] = inpolygon(real(radp{k}),imag(radp{k}),real(ver{j}),imag(ver{j}));
            radp{k}(inm)=NaN+i*NaN;  radp{k}(onm) =NaN+i*NaN;
        end
        radpvc{k}  = radp{k}(abs(radp{k})>=0);
        iradpvc{k} = evalu(f,radpvc{k},'d');
        iradp{k}  = NaN(size(radp{k}))+i*NaN(size(radp{k}));
        iradp{k}(abs(radp{k})>=0)  =iradpvc{k};        
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
for k=1:n_of_cr
    plot(real(cirp{k}),imag(cirp{k}),'-r');
end
for k=1:n_of_rd
    plot(real(radp{k}),imag(radp{k}),'-b');
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
for k=1:n_of_cr
    plot(real(icirp{k}),imag(icirp{k}),'-r');
end
for k=1:n_of_rd
    plot(real(iradp{k}),imag(iradp{k}),'-b');
end
% ax = gca;
% ax.XAxisLocation = 'origin';
% ax.YAxisLocation = 'origin';
%%
end