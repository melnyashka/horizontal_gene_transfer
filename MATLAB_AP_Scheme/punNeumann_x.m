%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fct de decalage vers la droite %%
%%        v(i)=u(i+1)             %%
%%        v(imax)=cl              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function v=punNeumann_x(u)
% 
% imax=size(u,2);
% cl=u(:,imax);
% v=0*u;
% v(:,1:imax-1)=u(:,2:imax);
% v(:,imax)=cl;

%le d�calage vers la droite devient un d�calage vers le bas
%Conditions de neumann
function v=punNeumann_x(u)

imax=size(u,1);
%lgn=u(imax,:);
lgn=2*u(imax,:)-u(imax-1,:);
%lgn=0*u(imax,:)+1;
v=0*u;
v(imax,:)=lgn;
v(1:imax-1,:)=u(2:imax,:);