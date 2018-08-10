function t=tau(X,Y)

%Fontion tau qui intervient dans l'EDP. Plusieurs possibilités à regarder.
%Taux à choisir : 
taux=2.2;   

%Pour les Heavyside :
%>= 0.5, extinction, \in[0.1, ] : oscillations

%Pour arctan :  
%0.1 : population qui stagne à ep=10^(-6)
%0.5 :  Si il y a oscillations à ep=10^(-6), elles sont lentes. En T=1.5 ce
%n'est pas clair
%0.6 : idem
%0.9 : idem
%1.5 : en t=15, c'est clair qu'il n'y a pas d'oscillations et pas
%d'extinction
%3 : 


%possibilité 1 : 
%x et y sont deux vecteurs ou matrices qui n'ont pas forcément la même
%taille. 
%t_{i,j}=tau(x_i,y_j)


%Possibilité 2 : 
%x et y sont deux vecteurs ou matrices qui ont la même taille. 
%t_{i,j} = tau(x_{i,j},y_{i,j}

% a priori, on va utiliser la possibilité 2:
%Avec la fonction indicatrice par exemple


%IMPORsTANT THINGS BELOW

%t=taux.*(X>Y)-taux.*(Y>X);
t=taux*tanh(X-Y);
