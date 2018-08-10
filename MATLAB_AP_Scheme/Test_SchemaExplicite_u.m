clear;
clf();

%Test du schéma explicite pour u
xmin=-2;
xmax=2;
Nx=200;

Tmax=1;
Nt=10000;

ep=10^(-5);


naissance='birth';
mort='death';
noyau='m';
taux='tau';
u0='u_init';

[t,dt,x,dx,u,U,f,F,rho,RHO]=feval('SchemaExplicite_u',ep,Tmax,Nt,xmin,xmax,Nx,u0,naissance,mort,taux,noyau);


%Tracés de la solution du schéma explicite sans la comparer à rien
%{
%Tracé de u juste au temps final
figure(1)
plot(x,u,'+-',x,U(:,1),'--','Linewidth',2,'Markersize',8)
xlabel('\fontsize{26} x')
legend('\fontsize{26} u(T_{max},x)','\fontsize{26} u(0,x)','Location','Best')
title(['\fontsize{26} T_{max} = ' num2str(Tmax) ', \Delta x  = ' num2str(dx) ', \Delta t = ' num2str(dt)])

%Tracé de f juste au temps final
figure(2)
plot(x,f,'+-',x,F(:,1),'--','Linewidth',2,'Markersize',8)
xlabel('\fontsize{26} x')
legend('\fontsize{26} f(T_{max},x)','\fontsize{26} f(0,x)','Location','Best')
title(['\fontsize{26} T_{max} = ' num2str(Tmax) ', \Delta x  = ' num2str(dx) ', \Delta t = ' num2str(dt)])

%Tracé de rho en fonction du temps
figure(3)
plot(t,RHO,'-+','Linewidth',2,'Markersize',8)
xlabel('\fontsize{26} t')
ylabel('\fontsize{26} \rho(t)')
title(['\fontsize{26} T_{max} = ' num2str(Tmax) ', \Delta x  = ' num2str(dx) ', \Delta t = ' num2str(dt)])

%Video de u et f en fonction du temps
disp('Appuyer sur une touche pour lancer la video')


figure(4)
subplot(1,2,1)
plot(x,u,'+-',x,U(:,1),'--','Linewidth',2,'Markersize',8)
xlabel('\fontsize{26} x')
legend('\fontsize{26} u(T_{max},x)','\fontsize{26} u(0,x)','Location','Best')
title(['\fontsize{26} T_{max} = ' num2str(Tmax) ', \Delta x  = ' num2str(dx) ', \Delta t = ' num2str(dt)])

subplot(1,2,2)
plot(x,f,'+-',x,F(:,1),'--','Linewidth',2,'Markersize',8)
xlabel('\fontsize{26} x')
legend('\fontsize{26} f(T_{max},x)','\fontsize{26} f(0,x)','Location','Best')
title(['\fontsize{26} T_{max} = ' num2str(Tmax) ', \Delta x  = ' num2str(dx) ', \Delta t = ' num2str(dt)])
pause

for i=1:50:length(t)
    
    subplot(1,2,1)
    plot(x,U(:,i),'+-',x,U(:,1),'--','Linewidth',2,'Markersize',8)
    xlabel('\fontsize{26} x')
    legend('\fontsize{26} u(t,x)','\fontsize{26} u(0,x)','Location','Best')
    title(['\fontsize{26} t = ' num2str(t(i)) ', \Delta x  = ' num2str(dx)])
    
    subplot(1,2,2)
    plot(x,F(:,i),'+-',x,F(:,1),'--','Linewidth',2,'Markersize',8)
    xlabel('\fontsize{26} x')
    legend('\fontsize{26} f(t,x)','\fontsize{26} f(0,x)','Location','Best')
    title(['\fontsize{26}  \Delta t = ' num2str(dt)])
    pause(0.01)
    
end
%}

%Tracés pour ep petit en comparant au schéma limite pour u
[tlim,dtlim,xlim,dxlim,ulim,Ulim,vlim,grandVlim,Xtlim]=feval('SchemaLimite_u',Tmax,Nt,xmin,xmax,Nx,u0,naissance,mort,taux,noyau);

[minU,places]=min(U);
Xt=x(places);

figure(1)
%subplot(1,2,1)
plot(x,u,'+-',x,ulim,'--',x,U(:,1),'o','Linewidth',2,'Markersize',8)
xlabel('\fontsize{26} x')
ylabel('\fontsize{26} u')
legend('\fontsize{26} u_{exp}(T_{max},x)','\fontsize{26} u_{lim}(T_{max},x)','\fontsize{26} u(0,x)','Location','Best')
title(['\fontsize{26}T_{max} = ' num2str(Tmax) ', \Delta x = ' num2str(dx) ])


figure(2)
subplot(1,3,1)
imagesc([0,Tmax],[xmin,xmax],Ulim)
xlabel('\fontsize{26} t')
ylabel('\fontsize{26} x')
title('\fontsize{26} u_{lim}(t,x)')

subplot(1,3,2)
imagesc([0,Tmax],[xmin,xmax],U)
xlabel('\fontsize{26} t')
ylabel('\fontsize{26} x')
title('\fontsize{26} u_{exp}(t,x)')

subplot(1,3,3)
MaxU=ones(length(x),1)*max(abs(Ulim));
imagesc([0,Tmax],[xmin,xmax],abs(U-Ulim)./MaxU)
xlabel('\fontsize{26} t')
ylabel('\fontsize{26} x')
title('\fontsize{26} |u_{exp}(t,x)-u_{lim}(t,x)|/|u_{lim}(t,.)|_{\infty}')


figure(3)
%subplot(2,2,1)
plot(t,Xt,'+-',t,Xtlim,'--','Linewidth',2,'MarkerSize',8)
xlabel('\fontsize{26} t')
ylabel('\fontsize{26} X(t) ')
legend('\fontsize{26} X_{exp}(t)','\fontsize{26} X_{lim}(t)')
title(['\fontsize{26} \Delta x = ' num2str(dx) ', \Delta t = ' num2str(dt)])





figure(4)
%subplot(1,2,1)
plot(x,u,'+-',x,ulim,'--',x,U(:,1),'o','Linewidth',2,'Markersize',8)
xlabel('\fontsize{26} x')
legend('\fontsize{26} u_{exp}(T_{max},x)','\fontsize{26} u_{lim}(T_{max},x)','\fontsize{26} u(0,x)','Location','Best')
title(['\fontsize{26} T_{max} = ' num2str(Tmax) ', \Delta x  = ' num2str(dx) ', \Delta t = ' num2str(dt)])


disp('Appuyer sur une touche pour afficher la video')
pause

for i=1:50:length(t)
    
    %subplot(1,2,1)
    plot(x,U(:,i),'+-',x,Ulim(:,i),'--',x,U(:,1),'o','Linewidth',2,'Markersize',8)
    xlabel('\fontsize{26} x')
    legend('\fontsize{26} u_{exp}(t,x)','\fontsize{26} u_{lim}(t,x)','\fontsize{26} u(0,x)','Location','Best')
    title(['\fontsize{26} t = ' num2str(t(i)) ', \Delta x  = ' num2str(dx)])
    
    pause(0.01)
    
end






