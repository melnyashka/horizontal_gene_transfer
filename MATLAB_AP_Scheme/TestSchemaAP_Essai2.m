clear;
clf();

%Script test du SchemaAP_Essai1 : le schéma tout seul, comparé au schéma
%pour f et comparé au schéma limite
ep=1;


xmin=-3;
xmax=3;
Nx=300;

Tmax=15;  %1.2;
Nt=30000; %12000;



naissance='birth';
mort='death';
noyau='m';
taux='tau';
u0='u_init';


[t,dt,x,dx,u,U,rho,RHO,f,F]=feval('SchemaAP_Essai2',ep,Tmax,Nt,xmin,xmax,Nx,u0,naissance,mort,taux);


%Tracés avec SchemaAP_Essai1 tout seul :
%{
figure(1)
subplot(2,2,1)
plot(x,u,'+-',x,U(:,1),'--','Linewidth',2,'Markersize',8)
xlabel('\fontsize{26} x')
ylabel('\fontsize{26} u')
legend('\fontsize{26} u(T_{max},x)','\fontsize{26} u(0,x)','Location','Best')
title(['\fontsize{26}T_{max} = ' num2str(Tmax) ', \Delta x = ' num2str(dx) ])

subplot(2,2,2)
plot(x,f,'+-',x,F(:,1),'--','Linewidth',2,'Markersize',8)
xlabel('\fontsize{26} x')
legend('\fontsize{26} f(T_{max},x)','\fontsize{26} f(0,x)','Location','Best')
title(['\fontsize{26}  \Delta t = ' num2str(dt)])

subplot(2,2,3)
imagesc([0,Tmax],[xmin,xmax],U)
xlabel('\fontsize{26} t')
ylabel('\fontsize{26} x')
title('\fontsize{26} u(t,x)')


subplot(2,2,4)
imagesc([0,Tmax],[xmin,xmax],F)
xlabel('\fontsize{26} t')
ylabel('\fontsize{26} x')
title('\fontsize{26} f(t,x)')


figure(2)
%subplot(2,2,1)
plot(t,RHO,'+-','Linewidth',2,'MarkerSize',8)
xlabel('\fontsize{26} t')
ylabel('\fontsize{26} \rho(t) ')
title(['\fontsize{26} \Delta x = ' num2str(dx) ', \Delta t = ' num2str(dt)])





figure(3)
subplot(1,2,1)
plot(x,u,'+-',x,U(:,1),'--','Linewidth',2,'Markersize',8)
xlabel('\fontsize{26} x')
legend('\fontsize{26} u(T_{max},x)','\fontsize{26} u(0,x)','Location','Best')
title(['\fontsize{26} T_{max} = ' num2str(Tmax) ', \Delta x  = ' num2str(dx) ', \Delta t = ' num2str(dt)])

subplot(1,2,2)
plot(x,f,'+-',x,F(:,1),'--','Linewidth',2,'Markersize',8)
xlabel('\fontsize{26} x')
legend('\fontsize{26} f(T_{max},x)','\fontsize{26} f(0,x)','Location','Best')
%title(['\fontsize{26} T_{max} = ' num2str(Tmax) ', \Delta x  = ' num2str(dx) ', \Delta t = ' num2str(dt)])

disp('Appuyer sur une touche pour afficher la video')
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

%Tracés pour ep grand en comparant au schéma explicite pour u
%{
[texp,dtexp,xexp,dxexp,uexp,Uexp,fexp,Fexp,rhoexp,RHOexp]=feval('SchemaExplicite_u',ep,Tmax,Nt,xmin,xmax,Nx,u0,naissance,mort,taux,noyau);

figure(1)
subplot(1,2,1)
plot(x,u,'+-',x,uexp,'--',x,U(:,1),'o','Linewidth',2,'Markersize',8)
xlabel('\fontsize{26} x')
ylabel('\fontsize{26} u')
legend('\fontsize{26} u_{AP}(T_{max},x)','\fontsize{26} u_{exp}(T_{max},x)','\fontsize{26} u(0,x)','Location','Best')
title(['\fontsize{26}T_{max} = ' num2str(Tmax) ', \Delta x = ' num2str(dx) ])

subplot(1,2,2)
plot(x,f,'+-',x,fexp,'--',x,F(:,1),'o','Linewidth',2,'Markersize',8)
xlabel('\fontsize{26} x')
legend('\fontsize{26} f_{AP}(T_{max},x)','\fontsize{26} f_{exp}(T_{max},x)','\fontsize{26} f(0,x)','Location','Best')
title(['\fontsize{26}  \Delta t = ' num2str(dt)])

figure(2)
subplot(2,3,1)
imagesc([0,Tmax],[xmin,xmax],Uexp)
xlabel('\fontsize{26} t')
ylabel('\fontsize{26} x')
title('\fontsize{26} u_{exp}(t,x)')

subplot(2,3,2)
imagesc([0,Tmax],[xmin,xmax],U)
xlabel('\fontsize{26} t')
ylabel('\fontsize{26} x')
title('\fontsize{26} u_{AP}(t,x)')

subplot(2,3,3)
MaxU=ones(length(x),1)*max(abs(Uexp));
imagesc([0,Tmax],[xmin,xmax],abs(U-Uexp)./MaxU)
xlabel('\fontsize{26} t')
ylabel('\fontsize{26} x')
title('\fontsize{26} |u_{AP}(t,x)-u_{exp}(t,x)|/|u_{exp}(t,.)|_{\infty}')

subplot(2,3,4)
imagesc([0,Tmax],[xmin,xmax],Fexp)
xlabel('\fontsize{26} t')
ylabel('\fontsize{26} x')
title('\fontsize{26} f_{exp}(t,x)')

subplot(2,3,5)
imagesc([0,Tmax],[xmin,xmax],F)
xlabel('\fontsize{26} t')
ylabel('\fontsize{26} x')
title('\fontsize{26} f_{AP}(t,x)')

subplot(2,3,6)
MaxF=ones(length(x),1)*max(abs(Fexp));
imagesc([0,Tmax],[xmin,xmax],abs(F-Fexp)./MaxF)
xlabel('\fontsize{26} t')
ylabel('\fontsize{26} x')
title('\fontsize{26}| f_{AP}(t,x)-f_{exp}(t,x)|/|u_{exp}(t,.)|_{\infty}')



figure(3)
%subplot(2,2,1)
plot(t,RHO,'+-',t,RHOexp,'--','Linewidth',2,'MarkerSize',8)
xlabel('\fontsize{26} t')
ylabel('\fontsize{26} \rho(t) ')
legend('\fontsize{26} \rho_{AP}(t)','\fontsize{26} \rho_{exp}(t)')
title(['\fontsize{26} \Delta x = ' num2str(dx) ', \Delta t = ' num2str(dt)])





figure(4)
subplot(1,2,1)
plot(x,u,'+-',x,uexp,'--',x,U(:,1),'o','Linewidth',2,'Markersize',8)
xlabel('\fontsize{26} x')
legend('\fontsize{26} u_{AP}(T_{max},x)','\fontsize{26} u_{exp}(T_{max},x)','\fontsize{26} u(0,x)','Location','Best')
title(['\fontsize{26} T_{max} = ' num2str(Tmax) ', \Delta x  = ' num2str(dx) ', \Delta t = ' num2str(dt)])

subplot(1,2,2)
plot(x,f,'+-',x,f,'--',x,F(:,1),'o','Linewidth',2,'Markersize',8)
xlabel('\fontsize{26} x')
legend('\fontsize{26} f_{AP}(T_{max},x)','\fontsize{26} f_{exp}(T_{max},x)','\fontsize{26} f(0,x)','Location','Best')
%title(['\fontsize{26} T_{max} = ' num2str(Tmax) ', \Delta x  = ' num2str(dx) ', \Delta t = ' num2str(dt)])

disp('Appuyer sur une touche pour afficher la video')
pause

for i=1:50:length(t)
    
    subplot(1,2,1)
    plot(x,U(:,i),'+-',x,Uexp(:,i),'--',x,U(:,1),'o','Linewidth',2,'Markersize',8)
    xlabel('\fontsize{26} x')
    legend('\fontsize{26} u_{AP}(t,x)','\fontsize{26} u_{exp}(t,x)','\fontsize{26} u(0,x)','Location','Best')
    title(['\fontsize{26} t = ' num2str(t(i)) ', \Delta x  = ' num2str(dx)])
    
    subplot(1,2,2)
    plot(x,F(:,i),'+-',x,Fexp(:,i),'--',x,F(:,1),'o','Linewidth',2,'Markersize',8)
    xlabel('\fontsize{26} x')
    legend('\fontsize{26} f_{AP}(t,x)','\fontsize{26} f_{exp}(t,x)','\fontsize{26} f(0,x)','Location','Best')
    title(['\fontsize{26}  \Delta t = ' num2str(dt)])
    pause(0.01)
    
end
%}

%Tracés pour ep petit en comparant au schéma limite pour u
%[tlim,dtlim,xlim,dxlim,ulim,Ulim,vlim,grandVlim,Xtlim]=feval('SchemaLimite_u',Tmax,Nt,xmin,xmax,Nx,u0,naissance,mort,taux,noyau);
%grandrholim=0*tlim;
%rholim=0;
[tlim,dtlim,xlim,dxlim,ulim,Ulim,Xtlim,rholim,grandrholim]=feval('SchemaLimite_Implicite_u',Tmax,Nt,xmin,xmax,Nx,u0,naissance,mort,taux,noyau);

[minU,places]=min(U);
Xt=x(places);

figure(1)
subplot(2,1,1)
plot(x,U(:,length(t)),'+-',x,Ulim(:,length(t)),'--',x,U(:,1),'o','Linewidth',2,'Markersize',8)
xlabel('\fontsize{26} x')
ylabel('\fontsize{26} u')
legend('\fontsize{26} u_{AP}(T_{max},x)','\fontsize{26} u_{lim}(T_{max},x)','\fontsize{26} u(0,x)','Location','Best')
title(['\fontsize{26}T_{max} = ' num2str(Tmax) ', \Delta x = ' num2str(dx) ', \epsilon = ' num2str(ep) ])

subplot(2,1,2)
plot(t,RHO,'+-',t,grandrholim,'--','Linewidth',2,'Markersize',8)
xlabel('\fontsize{26} t')
ylabel('\fontsize{26} \rho(t)')
legend('\fontsize{26} \rho_{AP}','\fontsize{26} \rho_{lim}')

figure(2)
subplot(1,3,1)
imagesc([0,Tmax],[xmin,xmax],Ulim)
xlabel('\fontsize{26} t')
ylabel('\fontsize{26} x')
title('\fontsize{26} u_{lim}(t,x)')

subplot(1,3,2)
imagesc([0,Tmax],[xmin,xmax],real(U))
xlabel('\fontsize{26} t')
ylabel('\fontsize{26} x')
title('\fontsize{26} u_{AP}(t,x)')

subplot(1,3,3)
MaxU=ones(length(x),1)*max(abs(Ulim));
imagesc([0,Tmax],[xmin,xmax],abs(U-Ulim)./MaxU)
xlabel('\fontsize{26} t')
ylabel('\fontsize{26} x')
title('\fontsize{26} |u_{AP}(t,x)-u_{lim}(t,x)|/|u_{lim}(t,.)|_{\infty}')


figure(3)
%subplot(2,2,1)
plot(t,Xt,'+-',t,Xtlim,'--','Linewidth',2,'MarkerSize',8)
xlabel('\fontsize{26} t')
ylabel('\fontsize{26} X(t) ')
legend('\fontsize{26} X_{AP}(t)','\fontsize{26} X_{lim}(t)')
title(['\fontsize{26} \Delta x = ' num2str(dx) ', \Delta t = ' num2str(dt) ', \epsilon = ' num2str(ep)])





figure(4)
%subplot(1,2,1)
plot(x,real(U(:,length(t))),'+-',x,Ulim(:,length(t)),'--',x,real(U(:,1)),'o','Linewidth',2,'Markersize',8)
xlabel('\fontsize{26} x')
legend('\fontsize{26} u_{AP}(T_{max},x)','\fontsize{26} u_{lim}(T_{max},x)','\fontsize{26} u(0,x)','Location','Best')
title(['\fontsize{26} T_{max} = ' num2str(Tmax) ', \Delta x  = ' num2str(dx) ', \Delta t = ' num2str(dt) ', \epsilon = ' num2str(ep)])


disp('Appuyer sur une touche pour afficher la video')
pause

for i=1:50:length(t)
    
    %subplot(1,2,1)
    plot(x,real(U(:,i)),'+-',x,Ulim(:,i),'--',x,U(:,1),'o','Linewidth',2,'Markersize',8)
    xlabel('\fontsize{26} x')
    legend('\fontsize{26} u_{AP}(t,x)','\fontsize{26} u_{lim}(t,x)','\fontsize{26} u(0,x)','Location','Best')
    title(['\fontsize{26} t = ' num2str(t(i)) ', \Delta x  = ' num2str(dx) ', \epsilon = ' num2str(ep)])
    
    pause(0.01)
    
end






















