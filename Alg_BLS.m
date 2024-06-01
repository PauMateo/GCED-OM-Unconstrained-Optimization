


f = @(x) x(1)*exp(-x(1)^2-x(2)^2);
g = @(x) exp(-x(1)^2-x(2)^2)*[(1 - 2*x(1)^2); (-2*x(1)*x(2))];x = [-1;1]; d = [3;-3];
almax = 1.0; almin = 10^-6; rho = 0.5; c1 = 0.1; c2 = 0.5; iW = 2;
rho = 0.99;

[al, iWout] = uo_BLS(x, d, f, g, almax, almin, rho, c1, c2, iW)

% [start] Alg. BLS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% [end] Alg. BLS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








