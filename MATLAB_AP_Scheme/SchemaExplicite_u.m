function [t,dt,x,dx,u,U,f,F,rho,RHO]=SchemaExplicite_u(ep,Tmax,Nt,xmin,xmax,Nx,u0,naissance,mort,taux,noyau)


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
%y=x';
%dy=dx;
ymax=8;
Ny=ceil(2*ymax/dx)+1;
dy=2*ymax/Ny;
y=-ymax:dy:ymax;



X=x*ones(1,length(y));
Y=ones(length(x),1)*y;

%to=feval(taux,X,Y);
to=feval(taux,x*ones(1,length(x)),ones(length(x),1)*x');
m=feval(noyau,y,dy);
mcarre=ones(length(x),1)*m;

%Moralité : pour intégrer, on somme les colonnes entre elles : sum(f,2)*dy
%(x_i,y_j)


for i=1:(length(t)-1)

    %Calcul de u 
    % Calcul de l'intégrale avec le b. Il va falloir se poser la question
    % des valeurs à mettre au bord. Pour l'instant on fixe : 
    % u(x) = u(xmin) si x<xmin 
    %u(x)= u(xmax) si x>xmax
    %{
    interpgrid=x'-ep*y;
    interpgrid=interpgrid.*(interpgrid>=xmin).*(interpgrid<=xmax)+xmin.*(interpgrid<xmin)+xmax.*(interpgrid>xmax);  %Pour gérer les conditions de bord.
    uinterp=interp1(x,u,interpgrid,'linear');  %Interpolation linéaire pour des questions de stabilité. Normalement sort en ligne.
    I1=(u*ones(1,length(y))-ones(length(x),1)*uinterp)/ep;
    b=feval(naissance,ones(length(x),1)*uinterp);
    I1=b.*exp(I1).*(ones(length(x),1)*m);
    I1=sum(I1,2)*dy;
    %}
    
    %Calcul de l'intégrale avec le b en mettant du Neumann au bord
    %{
    %{
    pentegauche=(u(2)-u(1))/dx;
    pentedroite=(u(length(x))-u(length(x)-1))/dx;
    interpgrid=x'-ep*y;  %Faux !!!
    interpgrid_bis=interpgrid.*(interpgrid>=xmin).*(interpgrid<=xmax)+xmin.*(interpgrid<xmin)+xmax.*(interpgrid>xmax);  %Pour gérer les conditions de bord.
    uinterp=interp1(x,u,interpgrid_bis,'linear');  %Interpolation linéaire pour des questions de stabilité. Normalement sort en ligne.
    uinterp=uinterp.*(interpgrid>=xmin).*(interpgrid<=xmax)...  %On ne change rien pour les points qui étaient à l'intérieur
        +(interpgrid<xmin).*(u(1)+pentegauche*(interpgrid-xmin))...  %On prolonge à gauche en utilisant pente gauche
        +(interpgrid>xmax).*(u(length(x))+pentedroite*(interpgrid-xmax));
    I1=(u*ones(1,length(y))-ones(length(x),1)*uinterp)/ep;
    b=feval(naissance,ones(length(x),1)*interpgrid);
    I1=b.*exp(I1).*(ones(length(x),1)*m);
    I1=sum(I1,2)*dy;
    %}
    
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
    
    ucarre=u*ones(1,length(y));
    ucarrep1=punNeumann_x(ucarre);
    ucarrem1=munNeumann_x(ucarre);
    
    pentegauche=(u(2)-u(1))/dx;
    pentedroite=(u(length(x))-u(length(x)-1))/dx;
    
    %
    flux=(Y>0).*(test<=dx).*Y.*(ucarre-ucarrem1)/dy  + (Y<=0).*(test<=dx).*Y.*(ucarrep1-ucarre)/dy... %(3-4)
        +(test>dx).*(interpgrid>=xmin).*(interpgrid<=xmax).*((ucarre-uinterp)./(ep))... %(1)
        +(test>dx).*(interpgrid<xmin).*Y.*pentegauche  ...  %(2-cas 1)
        +(test>dx).*(interpgrid>xmax).*Y.*pentedroite    ;  %(2-cas 2)
    %
    
    %{
    flux=(uinterp>=xmin).*(uinterp<=xmax).*Y.*(ucarre-uinterp)/ep...
        +(uinterp<xmin).*(u(1)+pentegauche*(interpgrid-xmin))...
        +(uinterp>xmax).*(u(length(x))+pentedroite*(interpgrid-xmax));
    %}
    %disp('flux = ')
    %disp(flux)
    
    I1=bcarre*exp(flux).*mcarre;
    I1=sum(I1,2)*dy;
    %}
    
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
    
    pentegauche=(u(2)-u(1))/dx;
    pentedroite=(u(length(x))-u(length(x)-1))/dx;
    ugauche=pentegauche*(interpgrid-xmin)+u(1);
    udroite=pentedroite*(interpgrid-xmax)+u(length(x));
    
    
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
    
    
    I1=bcarre.*exp(flux).*mcarre;
    I1=sum(I1,2)*dy;  %Exactement le même calcul que dans le schema AP
    
    
   %calcul de l'intégrale avec le taux:
   %{
   I2=to.*(ones(length(x),1)*f')/rho;
   I2=sum(I2,2)*dx;
   %}
   minu=min(u);
    rhou=sum(exp(-(u-minu)/ep))*dx;
    fsurrhou=exp(-(u-minu)/ep)/rhou;
    
    I2=to.*(ones(length(x),1)*fsurrhou');
    I2=sum(I2,2)*dx;
   
   
   
   %Mise à jour de u :
   u=u+dt*d+dt*rho-dt*I1-dt*I2;
    
    %Mise à jour de f et rho :
    f=exp(-u/ep);
    rho=sum(f)*dx;
    
    U(:,i+1)=u;
    RHO(i+1)=rho;
    F(:,i+1)=f;
    
    
end

