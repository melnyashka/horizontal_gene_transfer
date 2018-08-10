function [t,dt,x,dx,u,U,Xt,rho,grandrho]=SchemaLimite_Implicite_u(Tmax,Nt,xmin,xmax,Nx,u0,naissance,mort,taux,noyau)


%Calcule la solution du problème limite à l'aide d'un schéma implicite pour
%l'équation de Hamilton-Jacobi sur u
%On calcule rho avec un min
%Si rho<0, on arrête les calculs et on dit qu'il y a eu extinction


[t,dt,x,dx]=feval('grids',Tmax,Nt,xmin,xmax,Nx);

Xt=zeros(1,length(t));
grandrho=zeros(1,length(t));
grandrho(1)=NaN;

%Initialisation
u=feval(u0,x);

[~,place]=min(u);
Xt(1)=x(place);

U=zeros(length(x),length(t));
U(:,1)=u;


b=feval(naissance,x);
d=feval(mort,x);


%Définition de la grille d'intégration : ce sera x aussi, mais en ligne
%y=x';
%dy=dx;

ymax=8;
Ny=ceil(2*ymax/dx)+1;
dy=2*ymax/Ny;
y=-ymax:dy:ymax;



%X=x*ones(1,length(y));
Y=ones(length(x),1)*y;
Yplus=Y.*(Y>0);
Ymoins=Y.*(Y<0);



m=exp(-y.^2/2);
cste=sum(m)*dy;
mcarre=ones(length(x),1)*m/cste;

%m=feval(noyau,x,dx);

bcarre=b*ones(1,length(y));
%mcarre=m*ones(1,length(y));

for i=1:(length(t)-1)
    
    if (grandrho(i)<=0)
        
        U(:,i+1)=0*u;
        grandrho(i+1)=0;
        
    else
        
        %Calcul de zgradv, noté flux. Avec un schéma upwind
        um1=feval('munNeumann_x',u);
        up1=feval('punNeumann_x',u);
        Um1=um1*ones(1,length(y));
        Up1=up1*ones(1,length(y));
        UU=u*ones(1,length(y));
        
        flux=Yplus.*(UU-Um1)/dx+Ymoins.*(Up1-UU)/dx;
        
        %Calcul de l'intégrale
        I1=bcarre.*exp(flux-(ones(length(x),1)*y).^2/2)/cste;%.*mcarre;
        I1=sum(I1,2)*dy;
        
        
        %Calcul de tau(x,X(t))=tauX
        tauX=feval(taux,x,Xt(i)*ones(length(x),1));
        
        %Calcul de rho :
        truc=u/dt+d-I1-tauX;
        rho=-min(truc);
        
        %Mise à jour de u
        u=u+dt*rho+dt*d-dt*I1-dt*tauX;
        
        [val,place]=min(u);
        Xt(i+1)=x(place);
        
        
        
        U(:,i+1)=u;
        grandrho(i+1)=rho;
        
    end
    
    
end

