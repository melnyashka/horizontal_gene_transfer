function u=u_init(x)

%u=(x).^2;

 u=(abs(x)<=1).*x.^2/2+(abs(x)>1).*(abs(x)-1/2);


%u=(abs(x)<=1).*x.^2+(abs(x)>1).*(2*abs(x)-1);

