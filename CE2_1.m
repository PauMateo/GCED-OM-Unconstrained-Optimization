% Problem
clear;clc;
f = @(x) x(1)^2 + x(2)^3 + 3*x(1)*x(2);
g = @(x) [ 2*x(1)+3*x(2) ; 3*x(2)^2 + 3*x(1)]; h = @(x) [];


% Input parameters.
x1 = [-3;-1];                                           % starting solution.
epsG= 10^-6; kmax= 1500;                                % Stopping criterion:
almax= 2; almin= 10^-3; rho=0.5;c1=0.01;c2=0.45; iW= 2; % Linesearch:
isd= 2;                                                 % Search direction.
icg= 1; irc= 2; nu = 0.1;                              % Parameters for CGM method.
delta = 0;                                              % Parameter for MNN-CMI method.

n=10;
[xk,dk,alk,iWk,betak,Hk,tauk] = uo_solve(x1,f,g,h,epsG,kmax,almax,almin,rho,c1,c2,iW,isd,icg,n,irc,nu,delta);

xo = []; xylim = [0,0,0,0]; logfreq = 1;
[gk,la1k,kappak,rk,Mk] = uo_solve_log(x1,f,g,h,epsG,kmax,almax,almin,rho,c1,c2,iW,isd,icg,irc,nu,delta,xk,dk,alk,iWk,betak,Hk,tauk,xo,xylim,logfreq);






