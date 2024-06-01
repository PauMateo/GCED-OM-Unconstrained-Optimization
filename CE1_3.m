

f   = @(x) x^2 + 2*sin(x);
df  = @(x) 2*x + 2*cos(x);
ddf = @(x) 2 - 2*sin(x);

%parametres:
almax = 1.0; almin = 0.01;
rho = 0.75; c1 = 0.1; c2 = 0.5;
x0 = 5;
x_opt = -0.739085;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRADIENT METHOD
disp('GRADIENT METHOD');
k=1;
x=x0;
xant = x;
while norm(df(x)) >= 10^-6
    d = -df(x);
    
    [al, iWout] = uo_BLS(x, d, f, df, almax, almin, rho, c1, c2, 1);
    x = x + d*al;

    rr = norm(x - x_opt)/norm(xant - x_opt);
    M = norm(x - x_opt)/norm(xant - x_opt)^2;
    xant = x;
    fprintf('%3d %7.4f %7.4f %3.1e %6.4f\n', k, x, norm(df(x)), rr, M);
    k = k+1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEWTON METHOD
disp('');
disp('NEWTON METHOD');
k=1;
x=x0;
xant = x;
while norm(df(x)) >= 10^-6
    d = -df(x) / ddf(x);
    [al, iWout] = uo_BLS(x, d, f, df, almax, almin, rho, c1, c2, 1);
    x = x + d*al;
    rr = norm(x - x_opt)/norm(xant - x_opt);
    M = norm(x - x_opt)/norm(xant - x_opt)^2;
    xant = x;
    fprintf('%3d %7.4f %7.4f %3.1e %6.4f\n', k, x, norm(df(x)), rr, M);
    k = k+1;
end



function [al,iWout] = uo_BLS(x,d,f,g,almax,almin,rho,c1,c2,iW)
    %iW=1 -> WC
    %iW=2 -> SWC
    phi = @(al) f(x+al*d);
    dphi = @(al) g(x+al*d)'*d;
    WC1 = @(al) phi(al) <= phi(0) + c1*dphi(0)*al;
    WC2 = @(al) dphi(al) >= c2*dphi(0);
    WC = @(al) WC1(al) & WC2(al);
    SWC2 = @(al) abs(dphi(al)) <= abs(c2*dphi(0));
    SWC  = @(al) WC1(al) & SWC2(al);
    al = almax;
    while (al >= almin) 
        if (iW==1 && WC(al)) || (iW==2 && SWC(al))
            break
        end
        al = rho * al;
    end
    
    if SWC(al) && (iW==2)
        iWout = 3;
    elseif WC(al) && (iW==1)
        iWout = 2;
    elseif WC1(al) && (iW==1)
        iWout = 1;
    else
        iWout = 0;
    end
end