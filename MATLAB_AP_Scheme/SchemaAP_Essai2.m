function [t,dt,x,dx,u,U,rho,RHO,f,F]=SchemaAP_Essai2(ep,Tmax,Nt,xmin,xmax,Nx,u0,naissance,mort,taux)


%Schéma AP en prenant l'idée de Vincent : calculer d'abord rho en implicite
%et en déduire u


[t,dt,x,dx]=feval('grids',Tmax,Nt,xmin,xmax,Nx);

%Initialisation
u=feval(u0,x);
rho=sum(exp(-u/ep))*dx;
f=exp(-u/ep);


F=zeros(length(x),length(t));
RHO=zeros(1,length(t));
U=zeros(length(x),length(t));
U(:,1)=u;
RHO(1)=rho;
F(:,1)=f;

%b=feval(naissance,x);
d=feval(mort,x);


%Définition de la grille d'intégration : 
%Ce doit être une ligne
%Pour l'instant on prend x, mais il va falloir le changer
%y=x';
%dy=dx;

ymax=8;
Ny=ceil(2*ymax/dx)+1;
dy=2*ymax/Ny;
y=-ymax:dy:ymax;


X=x*ones(1,length(y));
Y=ones(length(x),1)*y;


to=feval(taux,x*ones(1,length(x)),ones(length(x),1)*x');  %L'intégrale qui contient tau sera caculée sur la même grille que u 

%m=exp(-x.^2/2);
%cste=sum(m)*dx;
m=exp(-y.^2/2);
cste=sum(m)*dy;
%mcarre=ones(length(x),1)*m/cste;


tol=10^(-13);
Nmax=10000;

for i=1:(length(t)-1)
    
    
    if (rho<=0)
        
        U(:,i+1)=zeros(length(x),1);
        RHO(i+1)=0;
        f=zeros(length(x),1);  %Ne sert pas pour les calculs, juste pour la sortie
        F(:,i+1)=f;  %idem
        
    else
    
        
        minu=min(u);
        %Calcul du second membre A :
        %Calcul de l'intégrale avec tau : T
        
        rhou=sum(exp(-(u-minu)/ep))*dx;
        fsurrhou=exp(-(u-minu)/ep)/rhou;
        
        T=to.*(ones(length(x),1)*fsurrhou');
        T=sum(T,2)*dx;
        
        
        
        %Calcul de l'intégrale avec b : H
        interpgrid=X-ep*Y;
        bcarre=feval(naissance,interpgrid);
        %On fait :
        %   - 1 Interpolation linéaire de Matlab là où interpgrid \in [xmin,xmax] et |X-interpgrid|>dx
        %   - 2 Extrapolation linéaire à la main là où interpgrid \notin [xmin,xmax] et |X-interpgrid|>dx
        %   - 3 Taux d'accroissement upwind là où interpgrid \in [xmin,xmax] et
        %       |X-interpgrid|<=dx
        %   - 4 Taux d'accroissement upwind avec extrapolation linéaire à la main
        %       là où interpgrid \notin [xmin,xmax] et |X-interpgrid|<=dx;
        
        test=ep*abs(Y);
        interpgridborne=interpgrid.*(interpgrid>=xmin).*(interpgrid<=xmax)+xmin*(interpgrid<xmin)+xmax*(interpgrid>xmax);
        uinterp=interp1(x,u,interpgridborne,'linear');
        
        %Ynonzero=Y+(Y==0);  %vaut Y partout où Y est non nul et 1 sinon
        
        ucarre=u*ones(1,length(y));
        ucarrep1=punNeumann_x(ucarre);
        ucarrem1=munNeumann_x(ucarre);
        
        pentegauche=(u(2)-u(1))/dx;   %left slope
        pentedroite=(u(length(x))-u(length(x)-1))/dx;  %right slope
        ugauche=pentegauche*(interpgrid-xmin)+u(1);  %u(xmin-dx)
        udroite=pentedroite*(interpgrid-xmax)+u(length(x));  %u(xmax+dx)
        
        
        flux=(test>=dx).*(...
            (interpgrid>=xmin).*(interpgrid<=xmax).*(ucarre-uinterp)/ep    ...
            + (interpgrid<xmin).*(ucarre-ugauche)/ep    ....
            +(interpgrid>xmax).*(ucarre-udroite)/ep ...
            )...
            +(test<dx).*(test>0).*(    ...
            (interpgrid>=xmin).*(interpgrid<=xmax).*(  ...
            (Y>0).*(ucarre-ucarrem1).*Y/dx ...
            +(Y<0).*(ucarrep1-ucarre).*Y/dx ...
            )...
            +(interpgrid<xmin).*pentegauche.*Y ...
            + (interpgrid>xmax).*pentedroite.*Y ...
            )...
            +(test==0)*0;
        
        
        
        
        %{
    disp(max(max(flux)))
    disp(min(min(flux)))
    pause
            %}
            
            %disp('flux = ')
            %disp(flux)
            
            
            H=bcarre.*exp(flux-(ones(length(x),1)*y).^2/2)/cste;%.*mcarre;
            H=sum(H,2)*dy;  
            
            
            
            %Méthode de Newton pour le calcul de rho : sur la formule
            %exponentielle. Ca n'est stable que jusque ep=10^(-4)
            %{
    A=u+dt*d-dt*H-dt*T;  %C'est une colonne
    C=dx*sum(exp(-A/ep));  %C'est un scalaire

    
    fonc=@(x)  x.*exp(dt*x/ep)-C;
    invdfonc=@(x) ep*exp(-dt*x/ep)./(ep+dt*x);
    
    
    err=10;
    compteur=1;
    
    r=rho;
    
    while (err>tol)&&(compteur<Nmax)
        compteur=compteur+1;
        
        newr=r-invdfonc(r).*fonc(r);
        newr=newr*(newr>=0); %Pour forcer la positivité. Ne devrait pas servir.
        err=abs(newr-r);
        r=newr;
        
    end
    
   
    rho=r;
            %}
            
            %Méthode de Newton pour le calcul de rho : sur la formule
            %logarithmique. Pour tenter d'être stable pour ep petit
            %
            
            A=u+dt*d-dt*H-dt*T;  %C'est une colonne
            minA=min(A);
            C=ep*log(dx)-minA+ep*log(sum(exp(-(A-minA)/ep)));  %C'est un scalaire
            fonc=@(x)  C-ep*log(x)-dt*x;
            invdfonc=@(x) -x./(ep+dt*x);
            
            
            err=10;
            compteur=1;
            
            r=rho;
            
            while (err>tol)&&(compteur<Nmax)
                compteur=compteur+1;
                
                newr=r-invdfonc(r).*fonc(r);
                newr=newr;%*(newr>=0); %Pour forcer la positivité. Ne devrait pas servir.
                err=abs(newr-r);
                r=newr;
                
            end
            
            
            rho=r;
            
            %Calcul de u :
            u=dt*rho+A;
            
            %Vérification : on doit avoir rho =<exp(-u/ep)>  : C'est vérifié
            
            
            U(:,i+1)=u;
            RHO(i+1)=rho;
            f=exp(-u/ep);  %Ne sert pas pour les calculs, juste pour la sortie
            F(:,i+1)=f;  %idem
            
            
    end
    
    
end





























