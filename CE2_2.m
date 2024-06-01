%    PRIMEA PART DEL QUESTIONARI 2
% clear;
% % Problem
% n=20; seed = 123456;
% [f,g,h,xo]=uo_sconvQF(n,seed);
% % Input parameters.
% x1 = rand(n,1);                                % first iterate.
% epsG= 10^-6; kmax= 1500;                       % Stopping criterium:
% almax= 0; almin= 0; rho=0; c1=0; c2=0; iW= 0;  % Linesearch:
% isd= 1; icg= 0; irc= 0 ; nu = 0; delta =0;     % Search direction
% % Optimization (tauk and delta are only useful with SDM).
% [xk,dk,alk,iWk,betak,Hk,tauk] = uo_solve(x1,f,g,h,epsG,kmax,almax,almin,rho,c1,c2,iW,isd,icg,irc,nu,delta);
% % Optimization log
% xylim = []; logfreq = 1;
% [gk,la1k,kappak,rk,Mk] = uo_solve_log(x1,f,g,h,epsG,kmax,almax,almin,rho, c1,c2,iW,isd,icg,irc,nu,delta,xk,dk,alk,iWk,betak,Hk,tauk,xo,xylim,logfreq);
% 
% clf; % clear figure window
% niter=size(xk,2); plot(diag(dk(:,1:niter-2)'* dk(:,2:niter-1)))


%   SEGONA PART DEL QÜESTIONARI 2


clear;
seed = 123456;
for n = 100:200:1100
    [f,g,h,xo]=uo_sconvQF(n,seed);
    % Input parameters.
    x1 = rand(n,1);                                % first iterate.
    epsG= 10^-6; kmax= 1500;                       % Stopping criterium.
    almax= 0; almin= 0; rho=0; c1=0; c2=0; iW= 0;  % Linesearch.
    isd= 2; icg= 1; irc= 0; nu = 0; delta =0;      % Search direction.
    % Optimization:
    [xk,dk,alk,iWk,betak,Hk,tauk] = uo_solve(x1,f,g,h,epsG,kmax, ...
        almax,almin,rho,c1,c2,iW,isd,icg,irc,nu,delta);
    fprintf( "n = %2i, niter = %2i\n", n,size(xk,2));
end




