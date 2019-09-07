function   [zet,zetp,cent,rad] = cirmapb(et,etp,alpha,nv,iprec,gmresrestart,gmrestol,gmresmaxit,koebetol,koebemaxit)
% cirmapb.m
% Nasser, September 5, 2019
% This function implement the numerical method presented in
% M.M.S. Nasser, Fast Computation of the Circular Map, Comput. Methods 
% Funct. Theory 15 (2015) 187-223
% to compute the conformal mapping w=f(z) from the bounded multiply 
% connected domain G onto the circular region D with the normaliziation
% f(alpha)=0 and f'(alpha)>0. Let et be the 
% parametrization of the boundary of the polygonal domain G and etp be 
% the first derivative of the parametrization et.
% This function compute: 
% zet: the parametrization of the boundary of the circular domain D; 
%      zet=f(et)
% zetp: the first derivative of the parametrization zet
% cent: a vector of length m where cent(k) is the center of the circle 
%       C_k=f(\Gamma_k), k=1,2,...,m, where cent(m)=0.
% Rad: a vector of length m where rad(k) is the radius of the circle 
%       C_k=f(\Gamma_k), k=1,2,...,m, where rad(m)=1.
%  
%
%
m     =  length(nv);
zet   =  et; 
zetp  =  etp;
for k=1:m-1
    cent(k) = getInteriorPoint(et(1+sum(nv(1:k-1)):sum(nv(1:k)),1), etp(1+sum(nv(1:k-1)):sum(nv(1:k)),1));
end
cent(m) = alpha;
Error =  1;
total_iteration = 0;
while (Error> koebetol)    
    total_iteration = total_iteration+1;
    zetold = zet;
    for k=1:m-1
        etm   = zet (1+sum(nv(1:k-1)):sum(nv(1:k)),1);
        etmp  = zetp(1+sum(nv(1:k-1)):sum(nv(1:k)),1);
        zm    = cent(k);
        [S,Sp,c]  =  cmu(etm,etmp,zm,nv(k),iprec,gmresrestart,gmrestol,gmresmaxit);
        wn    =  exp(-i*S);
        wnpep =  -i.*Sp.*wn;        
        % at the boundary k
        zet (1+sum(nv(1:k-1)):sum(nv(1:k)),1) = wn;
        zetp(1+sum(nv(1:k-1)):sum(nv(1:k)),1) = wnpep;
        cent(k) = 0;
        % for the other boundaries
        zet_other  = [zet(1:sum(nv(1:k-1)),1);zet(1+sum(nv(1:k)):sum(nv),1)]; 
        zet_otherp = [zetp(1:sum(nv(1:k-1)),1);zetp(1+sum(nv(1:k)):sum(nv),1)]; 
        Cau   = (zet_other-zm).*(fcaunu(etm,etmp,wn./(etm-zm),zet_other.',nv(k),c).');
        Cau2  = (fcauep(etm,etmp,wnpep,zet_other.',nv(k),c).');
        Caup  =  zet_otherp.*Cau2;
        zet (1:sum(nv(1:k-1)),1)       = Cau(1:sum(nv(1:k-1)),1);
        zetp(1:sum(nv(1:k-1)),1)       = Caup(1:sum(nv(1:k-1)),1);
        zet (1+sum(nv(1:k)):sum(nv),1) = Cau(1+sum(nv(1:k-1)):end,1);
        zetp(1+sum(nv(1:k)):sum(nv),1) = Caup(1+sum(nv(1:k-1)):end,1);
        for j=1:m
            if(j~=k)
                cent(j)  = (cent(j)-zm).*(fcaunu(etm,etmp,wn./(etm-zm),cent(j),length(etm),c).');
            end
        end        
        
    end   
    %external boundaries
    k = m;
    etm     = zet (1+sum(nv(1:k-1)):sum(nv(1:k)),1);
    etmp    = zetp(1+sum(nv(1:k-1)):sum(nv(1:k)),1);
    zm      = cent(k);
    [S,Sp]  =  cmb(etm,etmp,zm,length(etm),iprec,gmresrestart,gmrestol,gmresmaxit);
    wn  =  exp(i*S);
    wnpep =  i.*Sp.*wn;
    % at the boundary k
    zet (1+sum(nv(1:k-1)):sum(nv(1:k)),1) = wn;
    zetp(1+sum(nv(1:k-1)):sum(nv(1:k)),1) = wnpep;
    cent(k) = fcaunb(etm,etmp,wn,cent(k));
    % for the other boundaries
    zet_other  = zet (1:sum(nv(1:m-1)),1);
    zet_otherp = zetp(1:sum(nv(1:m-1)),1);
    Cau   = (fcaunb(etm,etmp,wn,zet_other.').');
    Cau2  = (fcauep(etm,etmp,wnpep,zet_other.').');
    Caup  =  zet_otherp.*Cau2;    
    zet (1:sum(nv(1:m-1)),1) = Cau;    
    zetp(1:sum(nv(1:m-1)),1) = Caup;  
    for j=1:m-1
        cent(j)  = fcaunb(etm,etmp,wn,cent(j));
    end
    
    for (j=1:m)
        rad(j,1)     = sum(abs(zet(1+sum(nv(1:j-1)):sum(nv(1:j)),1)-cent(j)))/nv(j);
    end
    if m>1
        Erroro   =  Error;
        Error    = norm(zet-zetold,inf)
        if (total_iteration>=10 & Error > Erroro)
            disp(['Koebe iterative method truncated after  ',...
                 num2str(total_iteration),'  iterations, with error  ',...
                 num2str(Error)])
            break;
        end
    end
    if m==1
        return;
    end
    if (total_iteration>=koebemaxit)
        disp(['Error of Koebe iterative method  ',...
               num2str(Error),...
              ' Convegence failed after maximum iterations allowed'])          
        break;
    end
end
end
%%