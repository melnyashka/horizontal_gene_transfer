%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% fct de decalage vers la gauche %%
%%%        v(i)=u(i-1)             %%
%%%        v(1)=cl                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function v=munNeumann_x(u)
% 
% imax=size(u,2);
% cl=u(:,1);
% v=0*u;
% v(:,1)=cl;
% v(:,2:imax)=u(:,1:imax-1);


%le d�calage vers la gauche devient un d�calage vers le haut
function v=munNeumann_x(u)

imax=size(u,1);
%lgn=u(1,:);
lgn=2*u(1,:)-u(2,:);
%lgn=0*u(1,:);
v=0*u;
v(2:imax,:)=u(1:imax-1,:);
v(1,:)=lgn;