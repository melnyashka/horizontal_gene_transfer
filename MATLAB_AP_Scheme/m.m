function y=m(x,dx)

%Calcule le noyau gaussien en s'assurant que sur l'intervalle il est de
%somme 1

y=exp(-x.^2/2);
c=sum(y)*dx;
y=y/c;