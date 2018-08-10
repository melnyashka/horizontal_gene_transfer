function [t,dt,x,dx,u,U,f,F,rho,RHO]=SchemaExplicite_f(ep,Tmax,Nt,xmin,xmax,Nx,u0,naissance,mort,taux)

[t,dt,x,dx]=feval('grids',Tmax,Nt,xmin,xmax,Nx);

%Initialisation
u=feval(u0,x);
f=exp(-u/ep);
rho=sum(f)*dx;

RHO=zeros(1,length(t));
F=zeros(length(x),length(t));
U=zeros(length(x),length(t));
F(:,1)=f;
U(:,1)=u;
RHO(1)=rho;

%b=feval(naissance,x);
d=feval(mort,x);


%Définition de la grille d'intégration : ce sera x aussi, mais en ligne
y=x';
dy=dx;

X=x*ones(1,length(y));
Y=ones(length(x),1)*y;

to=feval(taux,X,Y);
%m=feval(noyau,y,dy);  %ne servira pas, on le code directement dedans
b=feval(naissance,Y);  %b est directement un carré

mep=exp(-x.^2/(2*ep^2));
cste=sum(mep)*dx;
Mep=exp(-(X-Y).^2/(2*ep^2))/cste;


for i=1:(length(t)-1)
    
    %Calcul de l'intégrale avec le b, sans faire de changement de
    %variables dedans : I1
    %A voir si les bornes de la convolution sont les bonnes
    fcarre=ones(length(x),1)*f';
    I1=b.*fcarre.*Mep;
    I1=sum(I1,2)*dy;  %doit être une colonne
    
    
    %Calcul de l'intégrale avec le tau  : I2
    I2=to.*fcarre/rho;
    I2=sum(I2,2)*dy;  %Doit aussi être une colonne
    

    
    %Mise à jour de f :
    f=f-dt*(d+rho).*f/ep+dt*I1/ep+dt*f.*I2/ep;
    rho=sum(f)*dx;
    
    F(:,i+1)=f;
    u=-ep*log(f+eps);
    U(:,i+1)=u;
    RHO(i+1)=rho;
    
end

