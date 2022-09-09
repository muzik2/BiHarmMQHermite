function [] = Hermite2Dc()
    c = 0.0025;
    [x,y,n,ti,bi,li,ri,ii] = CreateModel();

    N = numel(x);
    bct = zeros(1,N);
    bcvD = zeros(1,N);
    bcvN = zeros(1,N);
    kB = [ti(:);bi(:);ri(:);li(:)];
    Ndb = numel(kB);

    bct(kB) = 1;
    bcvD(kB) = 0;
    bcvN(ti(:)) = 1;

    A = zeros(N+Ndb);
    b = zeros(N+Ndb,1);

    idb = 1;
    for i=1:N
        x0 = x(i);
        y0 = y(i);
        [P,dPdx,dPdy,dPdxx,dPdyy,dPdxy,LP,dLPdx,dLPdy,BP,dBPdx,dBPdy] = MQBase(x0,y0,x,y,c);
        if bct(i)==1
            A(i,:) = [P(:)',dPdx(kB)'.*n(kB,1)'+dPdy(kB)'.*n(kB,2)'];
            A(N+idb,:) = [dPdx(:)'*n(i,1)+dPdy(:)'*n(i,2), (dPdxx(kB)'.*n(kB,1)'+dPdxy(kB)'.*n(kB,2)')*n(i,1)+(dPdxy(kB)'.*n(kB,1)'+dPdyy(kB)'.*n(kB,2)')*n(i,2)];    
            b(i) = bcvD(i);
            b(N+idb) = bcvN(i);
            idb = idb+1;
        else
            A(i,:) = [BP(:)',dBPdx(kB)'.*n(kB,1)'+dBPdy(kB)'.*n(kB,2)'];
            b(i) = bcvD(i);
        end
    end
    
    alpha = A\b;

    bb = [min(x(:)),max(x(:)),min(y(:)),max(y(:))];

    [X,Y] = meshgrid(linspace(bb(1),bb(2),101),linspace(bb(3),bb(4),101));
    Z = X*0;
    NN = numel(Z);
    for i=1:NN
        x0 = X(i);
        y0 = Y(i);
        [P,dPdx,dPdy] = MQBase(x0,y0,x,y,c);
        base = [P(:)',dPdx(kB)'.*n(kB,1)'+dPdy(kB)'.*n(kB,2)'];
        val = base*alpha;
        Z(i) = val;
    end
    
    [v,u] = gradient(Z);
    surf(X,Y,Z);

    %contourf(X,Y,Z);
    %contour(X,Y,Z);
    hold on
    quiver(X(:),Y(:),u(:),-v(:),"b");
    hold off

    daspect([1,1,1]);
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
    nx = 101;
    ny = 101;
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
    
    %{
    ti = find(abs(y(:)-h)<tol);
    bi = find(abs(y(:))<tol);
    ri = find(abs(x(:)-w)<tol & y(:)-tol>0 & y(:)+tol<h);
    li = find(abs(x(:))<tol & y(:)-tol>0 & y(:)+tol<h);
    %}

    n(ti,1) = 0;
    n(ti,2) = 1;
    n(bi,1) = 0;
    n(bi,2) = -1;

    n(li,1) = -1;
    n(li,2) = 0;
    n(ri,1) = 1;
    n(ri,2) = 0;
    
    




    ii = find((y(:)-tol)>0 & (y(:)+tol)<h & (x(:)-tol)>0 & (x(:)+tol)<w);

    plot(x(ti),y(ti),"cs",x(bi),y(bi),"rs",x(li),y(li),"kd",x(ri),y(ri),"bo",x(ii),y(ii),"k+");
    hold on
    quiver(x(k),y(k),n(k,1),n(k,2),"b");
    hold off
end

function [P,dPdx,dPdy,dPdxx,dPdyy,dPdxy,LP,dLPdx,dLPdy,BP,dBPdx,dBPdy] = MQBase(x0,y0,x,y,c)
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

    BP = (-8*c^4+8*c^2*(r2+r2).^2)./(P.^7);
    dBPdx = -(3*dx.*(-24*c^4+12*c^2*(r2+r2).^2))./(P.^9);
    dBPdy = -(3*dy.*(-24*c^4+12*c^2*(r2+r2).^2))./(P.^9);
end

