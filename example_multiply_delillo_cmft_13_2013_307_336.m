clear all
%%
ver{1}=[ 2.31+1.72i ;  3.61+1.57i ;  2.70+0.81i              ];
ver{2}=[ 3.73-0.94i ;  3.81-1.68i ;  2.41-1.81i ;  2.64-0.79i];
ver{3}=[ 0.71-0.54i ; -0.71-0.64i ; -0.38+0.85i ;  0.56+0.64i];
alpha = inf;
%%
% f=plgcirmap(ver,alpha);
f=plgcirmap(ver,alpha,ver{end}(end));
%%
plotmap(f);
%%
plotmap(f,'d','rec',21,21);
%% 
plotmap(f,'d','plr',21,21);
%% 
plotmap(f,'v','rec',21,21);
%%
plotmap(f,'v','plr',21,21);
%%