function [] = Hermite1D()
    nx = 5;
    x = linspace(0,2*pi,nx);
    n = [-1;x(2:end-1)'*0;1]';
    bk = [1,numel(x)];
    y = MyFun(x);

    N = numel(x);
    Nb = 2;

    R0 = zeros(N);
    b0 = zeros(N,1);

    for i=1:N
        x0 = x(i);
        R = RBFBase(x0,x);
        F = MyFun(x0);
        R0(i,:) = R;
        b0(i) = F;
    end
    
    Rdb = zeros(N,Nb);
    Rdb1 = zeros(Nb,N);
    for i=1:N
        x0 = x(i);
        [~,dR] = RBFBase(x0,x(bk));
        Rdb(i,:) = dR.*n(bk);
    end
    
    Rc = zeros(Nb);
    bc = zeros(Nb,1);
    for i=1:Nb
        x0 = x(bk(i));
        [~,dF] = MyFun(x0);
        bc(i) = dF.*n(bk(i));
        [~,~,ddR] = RBFBase(x0,x(bk));
        Rc(i,:) = ddR.*n(bk).*n(bk(i));
    end

    R = [R0,Rdb;Rdb',Rc];
    b = [b0;bc];
    %R = R0;
    %b = b0;

    alpha = R\b;
    alpha0 = R0\b0;

    X = linspace(x(1),x(end),1000);
    Y = X*0;
    Y0 = X*0;
    Y1 = X*0;
    NN = numel(X);
    for i=1:NN
        x0 = X(i);
        [R] = RBFBase(x0,x);
        [~,dR] = RBFBase(x0,x(bk));
        G = [R,dR.*n(bk)];
        %G = R;
        val0 = R*alpha0;
        Y0(i) = val0;
        val = G*alpha;
        Y(i) = val;
        F = MyFun(x0);
        Y1(i) = F;
    end
    plot(X,Y,"r-",X,Y0,"k--",x,y,"b+",X,Y1,"b:");

    %{
    x = linspace(-1,2);
    [R,dR,ddR] = RBFBase(0,x);
    plot(x,R,"r-",x,dR,"b.-",x,ddR,"k--");
    %}
end

function [F,dF] = MyFun(x)
    F = sin(x);
    dF = cos(x);
end

function [R,dR,ddR] = RBFBase(x0,x)
    c = 4.1;

    dx = x-x0;
    r2 = dx.^2;
    r = sqrt(r2);
    sind = find(r<1e-8);
    
    R = sqrt(r2(:)'+c^2);
    dR = dx(:)'.*R;
    ddR = c^2./R.^3;
    
    %R(sind) = 0;
    %dR(sind) = 0;
    %ddR(sind) = 3;
end


%{
function [R,dR,ddR] = RBFBase(x0,x)
    dx = x-x0;
    r2 = dx.^2;
    r = sqrt(r2);
    sind = find(r<1e-8);
    R = r2(:)'.*log(r(:)');
    dR = dx(:)'.*(1+log(r2(:)'));
    ddR = 3 + log(r2(:)');
    
    R(sind) = 0;
    dR(sind) = 0;
    ddR(sind) = 3;
end


function [R,dR,ddR] = RBFBase(x0,x)
    dx = x-x0;
    r2 = dx.^2;
    r = sqrt(r2);
    sind = find(r<1e-8);
    R = r2(:)'.*r(:)';
    dR = 3*dx(:)'.*r(:)';
    ddR = 3*r2(:)'./r(:)'+3.*r(:)';
    %dR(sind) = 0;
    ddR(sind) = 0;
end
%}