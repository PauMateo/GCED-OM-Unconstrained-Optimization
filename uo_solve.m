% [start] Function [uo_solve] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xk,dk,alk,iWk,betak,Hk,tauk] = uo_solve(x,f,g,h,epsG,kmax,almax,almin,rho,c1,c2,iW,isd,icg,irc_n,irc,nu,delta)
%
% Unconstraied optimization algorithm with first and second derivative methods.
% c1    : (WC1) coefficient.
% c2    : (WC2) ans (SWC2) coefficient.
% iW    : 0 = exact LS; 1= WC; 2= SWC; 3=WC+SDC; 4= BLSNW32
% isd   : search direction: 1=GM; 2=CGM; 3=BFGS; 4=NM; 5=MNM-SD; 6=MNM-CMI.
% icg   : CGM variant: 1=FR; 2=PR+; 
% irc   : Restart for the CGM: 0= no restart; 1=(RC1); 2=(RC2).
% irc_n     : cada quantes iteracions es fa un restart pel cas irc=1
% nu    : parameter nu for the (RC2)
% delta : parameter delta in the MNM-CMI.
%
n=size(x,1);
xk = x; dk = []; alk = []; iWk = []; Hk=zeros(n); betak=[]; tauk=[];
I = eye(n,n);
H = I;
k = 1;
ldescent = true;

RC1 = @(k) mod(k,n)==0;
RC2 = @(x,x_ant) (abs(g(x)'*g(x_ant)) / norm(g(x))^2) >= nu;
RC = @(irc,x,x_ant,k) (irc==1 && RC1(k)) || (irc==2 && RC2(x,x_ant));
Q = h(zeros(n,1)); %??

while norm(g(x)) > epsG & k < kmax & (ldescent | isd == 4) %............... main loop

    if isd == 1 % GM
        d = -g(x);

    elseif isd == 2 % CGM
        if k == 1 | RC(irc,x,x_ant,k)
            beta = 0;
            d = -g(x);
        else
            if icg == 1 %FR
                beta = (g(x)'*g(x)) /norm(g(x_ant))^2;    
            else %PR+
                beta = max(0,(g(x)'*( g(x)-g(x_ant) ))/norm(g(x_ant))^2);

            end
                    d = -g(x) + beta*d_ant;
        end
        betak = [betak,beta];

    elseif isd == 3 % BFGS
        d = - H * g(x);

        Hk(:,:,k)=H;
    elseif isd == 4 % NM
        d = -inv(h(x)) * g(x);

    elseif isd == 5 % MNM-SD
        [Q, A] = eig(h(x));
        for i = 1:length(h(x))
            A(i,i) = max (A(i,i), delta);
        end

        B = Q * A * Q';
        d = -inv(B)*g(x);

        Hk(:,:,k)=B;      
    elseif isd == 6 % MNM-CMI
        laUB = norm(h(x), 'fro');
        it = 0; tau=0;
        [R,p] = chol( h(x) + tau*eye(n));
        while p ~= 0
            it = it + 1;
            tau = (1.01 - 1/(2^it)) * laUB;
            [R,p] = chol( h(x) + tau*eye(n));
            
        end
        B = R'*R;
        d = -inv(B) * g(x);

        Hk(:,:,k)=B;
        tauk = [tauk, tau];
    end
    ldescent = d'*g(x) < 0;


    % LS %%%%%%%%%%%%%%%%%%%%%
    if isd == 4         % unit step length for the Newton method.
        al = 1;
        iWout=4;
    elseif iW == 0      % Exact line search
        al = -(g(x)'*d)/(d'*Q*d);
        iWout=3;
    elseif iW <= 3      % BLS, constant reduction.
        [al, iWout] = uo_BLS(x,d,f,g,almax,almin,rho,c1,c2,iW);
    end


    % Update %%%%%%%%%%%%%%%%%%%
    x_ant = x; d_ant = d;
    x = x + al*d; k=k+1;

    % BFGS update:
    if isd == 3
        s = x - x_ant;
        y = g(x) - g(x_ant);
        p = 1/(y'*s);
        H = (I - p*s*y')*H*(I - p*y*s') + p*(s*s') ;
    end

    xk = [xk,x]; dk = [dk,d]; alk = [alk,al]; iWk = [iWk,iWout];
    
end %............................................................ main loop
iWk=[iWk,NaN];
end
% [end] Function [uo_solve] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

