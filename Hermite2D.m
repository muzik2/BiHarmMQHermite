function [] = Hermite2D()
    %[x,y] = meshgrid(linspace(-2,2));
    c = 0.05;
    %[P,dPdx,dPdy,dPdxx,dPdyy,dPdxy] = MQBase(0,0,x,y,c);
    [x,y,n,ti,bi,li,ri,ii] = CreateModel();

    N = numel(x);
    bct = zeros(1,N);
    bcv = zeros(1,N);
    kN = [ti(:);bi(:)];
    kD = [ri(:);li(:)];
    Ndb = numel(kN);

    bct(kD) = 1;
    bcv(kD) = 0;
    p = x(kD);
    bct(kN) = 2;
    bcv(kN) = 1;
    bcv(ii) = 0;

    A = zeros(N+Ndb);
    b = zeros(N+Ndb,1);

    idb = 1;
    for i=1:N
        x0 = x(i);
        y0 = y(i);
        [P,dPdx,dPdy,dPdxx,dPdyy,dPdxy,LP,dLPdx,dLPdy] = MQBase(x0,y0,x,y,c);
        if bct(i)==1
            A(i,:) = [P(:)',dPdx(kN)'.*n(kN,1)'+dPdy(kN)'.*n(kN,2)'];
            b(i) = bcv(i);
        elseif bct(i)==2
            % bebebebe
            A(i,:) = [LP(:)',dLPdx(kN)'.*n(kN,1)'+dLPdy(kN)'.*n(kN,2)'];
            A(N+idb,:) = [dPdx(:)'*n(i,1)+dPdy(:)'*n(i,2), (dPdxx(kN)'.*n(kN,1)'+dPdxy(kN)'.*n(kN,2)')*n(i,1)+(dPdxy(kN)'.*n(kN,1)'+dPdyy(kN)'.*n(kN,2)')*n(i,2)];    
            b(i) = 0;
            b(N+idb) = bcv(i);
            idb = idb+1;
        else
            A(i,:) = [LP(:)',dLPdx(kN)'.*n(kN,1)'+dLPdy(kN)'.*n(kN,2)'];
            b(i) = bcv(i);
        end
    end
    
    alpha = A\b;

    bb = [min(x(:)),max(x(:)),min(y(:)),max(y(:))];

    [X,Y] = meshgrid(linspace(bb(1),bb(2),51),linspace(bb(3),bb(4),51));
    Z = X*0;
    NN = numel(Z);
    for i=1:NN
        x0 = X(i);
        y0 = Y(i);
        [P,dPdx,dPdy] = MQBase(x0,y0,x,y,c);
        base = [P(:)',dPdx(kN)'.*n(kN,1)'+dPdy(kN)'.*n(kN,2)'];
        val = base*alpha;
        Z(i) = val;
    end

    surf(X,Y,Z);
    %{
    surf(x,y,P);
    xlabel("X");
    ylabel("Y");
    zlabel("Z");
    %}
end

function [x,y,n,ti,bi,li,ri,ii] = CreateModel()
    w = 1;
    h = 1;
    nx = 11;
    ny = 11;
    [x,y] = meshgrid(linspace(0,w,nx),linspace(0,w,ny));

    N = numel(x);
    n = zeros(N,2);
    
    k = boundary(x(:),y(:),0.9);
    k = k(1:end-1);

    Nk = numel(k);
    for i=1:Nk 
        ia = i-1;
        ib = i+1;
        if(i==1)
          ia = Nk;
        elseif(i==Nk)
          ib = 1;
        end
        tvec = [x(k(ib))-x(k(ia)),y(k(ib))-y(k(ia))];
        nvec = [tvec(2),-tvec(1)]/norm(tvec);
        n(k(i),:) = nvec;
    end
    tol = 1e-7;
    ti = find(abs(y(:)-h)<tol & x(:)-tol>0 & x(:)+tol<w);
    bi = find(abs(y(:))<tol & x(:)-tol>0 & x(:)+tol<w);
    ri = find(abs(x(:)-w)<tol);
    li = find(abs(x(:))<tol);
    ii = find((y(:)-tol)>0 & (y(:)+tol)<h & (x(:)-tol)>0 & (x(:)+tol)<w);

    plot(x(ti),y(ti),"cs",x(bi),y(bi),"rs",x(li),y(li),"kd",x(ri),y(ri),"bo",x(ii),y(ii),"k+");
    hold on
    quiver(x(k),y(k),n(k,1),n(k,2),"b");
    hold off
end

function [P,dPdx,dPdy,dPdxx,dPdyy,dPdxy,LP,dLPdx,dLPdy] = MQBase(x0,y0,x,y,c)
    arguments
        x0 (1,1) {mustBeReal, mustBeFinite} = 0;
        y0 (1,1) {mustBeReal, mustBeFinite} = 0;
        x (:,:) {mustBeReal, mustBeFinite} = [0,1,2];
        y (:,:) {mustBeReal, mustBeFinite} = [0,1,2];
        c (1,1) {mustBeReal, mustBeFinite} = 2.7;
    end

    dx = x-x0;
    dy = y-y0;

    r2 = dx.^2+dy.^2;

    P = sqrt(dx.^2+dy.^2+c^2);
    dPdx = dx./P; 
    dPdy = dy./P;
    dPdxx = (c.^2+dy.^2)./(P.^3);
    dPdyy = (c.^2+dx.^2)./(P.^3);
    dPdxy = -(dx.*dy)./(P.^3);

    LP = (2*c^2+r2)./(P.^3);
    dLPdx = -(dx.*(4*c^2+r2))./(P.^5);
    dLPdy = -(dy.*(4*c^2+r2))./(P.^5);
end