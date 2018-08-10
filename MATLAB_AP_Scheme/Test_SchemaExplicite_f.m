clear;
clf();

%Test du schéma explicite pour f
xmin=-2;
xmax=2;
Nx=200;

Tmax=0.1;
Nt=100;

ep=1;


naissance='birth';
mort='death';
%noyau='m';
taux='tau';
u0='u_init';

[t,dt,x,dx,u,U,f,F,rho,RHO]=feval('SchemaExplicite_f',ep,Tmax,Nt,xmin,xmax,Nx,u0,naissance,mort,taux);

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


figure(5)
%contourf(t,x,F)
%surf(t,x,F)
imagesc([0,Tmax],[xmin,xmax],F)
xlabel('\fontsize{26} t')
ylabel('\fontsize{26} x')
title('\fontsize{26} f(t,x)')

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

for i=1:20:length(t)
    
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

