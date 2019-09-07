function   plotmap(f,a,b,no1,no2)
%
%
if nargin==1
    plotbd(f);
elseif nargin==3 & a=='v' & b=='plr'
    no1 = 15; % number of concentric circles to be plotted in the domain
    n_of_rd = 20; % number of radials to be plotted in the domain
    plotplri(f,no1,n_of_rd)
elseif nargin==5 & a=='v' & b=='plr'
    plotplri(f,no1,no2)
elseif nargin==3 & a=='v' & b=='rec'
    no1 = 15; % number of concentric circles to be plotted in the domain
    n_of_rd = 20; % number of radials to be plotted in the domain
    plotreci(f,no1,n_of_rd)
elseif nargin==5 & a=='v' & b=='rec'
    plotreci(f,no1,no2)
elseif nargin==3 & a=='d' & b=='rec'
    no1 = 15; % number of horizontal lines to be plotted in the domain
    no2 = 20; % number of vertical lines to be plotted in the domain
    plotrecd(f,no1,no2);
elseif nargin==5 & a=='d' & b=='rec'
    plotrecd(f,no1,no2);
elseif nargin==3 & a=='d' & b=='plr'
    no1 = 15; % number of horizontal lines to be plotted in the domain
    no2 = 20; % number of vertical lines to be plotted in the domain
    plotplrd(f,no1,no2);
elseif nargin==5 & a=='d' & b=='plr'
    plotplrd(f,no1,no2);
end    
% 
end