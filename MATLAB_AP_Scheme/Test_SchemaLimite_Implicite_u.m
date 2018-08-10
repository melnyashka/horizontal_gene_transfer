clear;
clf();

%Test du schéma limite implicite sur u (et pas v). 
%Où il y a calcul de rho
%Pour tau=0, il doit renvoyer la même chose que l'autre schéma limite. 
%On espère qu'il soit meilleur pour tau\neq 0
xmin=-2; %-4;  %-2
xmax=2; %4;  %2
Nx=400;%200;  %100

Tmax=1;%0.1;  %2
Nt=1000;%10000;  %400



naissance='birth';
mort='death';
noyau='m';
taux='tau';
u0='u_init';

[t,dt,x,dx,u,U,Xt,rho,grandrho]=feval('SchemaLimite_Implicite_u',Tmax,Nt,xmin,xmax,Nx,u0,naissance,mort,taux,noyau);

%Tracé de u juste au temps final
figure(1)
plot(x,u,'+-',x,U(:,1),'--','Linewidth',2,'Markersize',8)
xlabel('\fontsize{26} x')
legend('\fontsize{26} u(T_{max},x)','\fontsize{26} u(0,x)','Location','Best')
title(['\fontsize{26} T_{max} = ' num2str(Tmax) ', \Delta x  = ' num2str(dx) ', \Delta t = ' num2str(dt)])


figure(3)
subplot(2,1,1)
plot(t,Xt,'+-','Linewidth',2,'MarkerSize',8)
xlabel('\fontsize{26} t')
ylabel('\fontsize{26} X(t)')

subplot(2,1,2)
plot(t,grandrho,'+-','Linewidth',2,'Markersize',8)
xlabel('\fontsize{26} t')
ylabel('\fontsize{26} \rho(t)')

%Video de u 
disp('Appuyer sur une touche pour lancer la video')


figure(2)
%subplot(1,2,1)
plot(x,u,'+-',x,U(:,1),'--','Linewidth',2,'Markersize',8)
xlabel('\fontsize{26} x')
legend('\fontsize{26} u(T_{max},x)','\fontsize{26} u(0,x)','Location','Best')
title(['\fontsize{26} T_{max} = ' num2str(Tmax) ', \Delta x  = ' num2str(dx) ', \Delta t = ' num2str(dt)])


pause

for i=1:50:length(t)
    
    figure(2)
    %subplot(1,2,1)
    plot(x,U(:,i),'+-',x,U(:,1),'--','Linewidth',2,'Markersize',8)
    xlabel('\fontsize{26} x')
    legend('\fontsize{26} u(t,x)','\fontsize{26} u(0,x)','Location','Best')
    title(['\fontsize{26} t = ' num2str(t(i)) ', \Delta x  = ' num2str(dx)])
    pause(0.01)
    
    
end