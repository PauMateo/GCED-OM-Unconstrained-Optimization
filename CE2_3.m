clear;
% Rosenbrock function
f = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
g = @(x) [ -400*(x(2)-x(1)^2)*x(1) - 2*(1-x(1));  200*(x(2)-x(1)^2) ];
h = @(x) [ -400*x(2) + 1200*x(1)^2+2  -400*x(1); -400*x(1) 200  ];
x1 = [-1.; 2.0];
% Stopping conditions:
epsG = 10^-8; kmax = 20000;
% Line search
almax = 1; almin = 10^-3; rho=0.9; c1=0.01; c2=0.45; iW = 2;
% Dummy parameter for SDM.
delta = 0;


%%%%%%%%%%%%%%%
clc;
isd = 3;
icg = 2;
irc = 2;
nu = 0.5;
%%%%%%%%%%%%%%%
% isd   : search direction: 1=GM; 2=CGM; 3=BFGS; 4=NM; 5=MNM-SD; 6=MNM-CMI.
% icg   : CGM variant: 1=FR; 2=PR+; 
% irc   : Restart for the CGM: 0= no restart; 1=(RC1); 2=(RC2).

irc_n=10;
[xk,dk,alk,iWk,betak,Hk,tauk] = uo_solve(x1,f,g,h,epsG,kmax,almax,almin,rho,c1,c2,iW,isd,icg,irc_n,irc,nu,delta);

% Iterations log
xo=[1.0;1.0]; xylim = []; logfreq = @(xk) ceil(size(xk,2)/100);


[gk,la1k,kappak,rk,Mk] = uo_solve_log(x1,f,g,h,epsG,kmax,almax,almin,rho,c1,c2,iW,isd,icg,irc,nu,delta,xk,dk,alk,iWk,betak,Hk,tauk,xo,xylim,logfreq(xk));


