function [t,dt,x,dx]=grids(Tmax,Nt,xmin,xmax,Nx)

dt=Tmax/Nt;
t=0:dt:Tmax;

dx=(xmax-xmin)/Nx;
x=xmin:dx:xmax;
x=x';

