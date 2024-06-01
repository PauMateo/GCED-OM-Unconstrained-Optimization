% Problem and first iterate
f  = @(x) x(1)^2+x(2)^3+3*x(1)*x(2);
g  = @(x) [ 2*x(1)+3*x(2); 3*x(2)^2+3*x(1)];
h  = @(x) [ 2 , 3; 3 , 6*x(2)];
x1 = [-2;0];
% Input parameters.

epsG = sqrt(eps); kmax = 100;                               % Stopping criterium
almax = 1.0; almin = 10^-6; rho=0.5;c1=0.01;c2=0.9; iW = 1; % Linesearch
isd = 5; delta = 0.1;
icg=0;
irc=0;
nu=0;

irc_n=10;
[xk,dk,alk,iWk,betak,Hk,tauk] = uo_solve(x1,f,g,h,epsG,kmax,almax,almin,rho,c1,c2,iW,isd,icg,irc_n,irc,nu,delta);

% Iterations log
xo=[1.0;1.0]; xylim = []; logfreq = @(xk) ceil(size(xk,2)/100);
[gk,la1k,kappak,rk,Mk] = uo_solve_log(x1,f,g,h,epsG,kmax,almax,almin,rho,c1,c2,iW,isd,icg,irc,nu,delta,xk,dk,alk,iWk,betak,Hk,tauk,xo,xylim,logfreq(xk));

