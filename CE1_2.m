

clear; n=1000; seed = 1234; [f,g,h,xo,b]=uo_sconvQF(n,seed); Q=h(xo);
x  = zeros(n,1);
almax = 1.0;
almin = 0.01;
rho = 0.5;
c1 = 0.1;
c2 = 0.5;

warning off;
disp("ELS")
tic
[xresELS, k] = GUOA(x, f, g, Q, b, almax, almin, rho, c1, c2, 3);
toc
disp(k);

disp("")
disp("BLS-SWC")
tic
[xresBLS, k] = GUOA(x, f, g, Q, b, almax, almin, rho, c1, c2, 2);
toc
disp(k);



%%%%%%%%%%%%%%%%%%%%%%%%%%% random function generator %%%%%%%%%%%%%%%%%
function [f,g,h,xo,b] = uo_sconvQF(n,seed)
    rng(seed); Q=rand(n);
    [V,D] = eig(triu(Q)+triu(Q)');
    Q=V*diag(diag(max(D,1)))*V';
    b=rand(n,1); xo = Q^-1*b;
    f = @(x) x'*Q*x/2-b'*x; g = @(x) Q*x-b;
    h = @(x) Q;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, k] = GUOA(x, f, g, Q, b, almax, almin, rho, c1, c2, wi)
    % wi = 1 -> WC
    % wi = 2 -> SWC
    % wi = 3 -> exact line search
    k = 0;
    while norm(g(x)) >= 10^-6
        d = -g(x);
        if wi<3
            [al, iWout] = uo_BLS(x,d,f,g,almax,almin,rho,c1,c2,wi);
        else
            al = -(g(x)'*d)/(d'*Q*d);
        end
        x = x + al*d;
        k = k+1;

    end
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












