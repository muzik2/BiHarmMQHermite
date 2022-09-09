function [] = HermiteInterpolation()
   debug = ~true;
   na = 11;
   [xa,ya] = meshgrid(linspace(0,1,na));
   [k,n] = GetBoundaryAndNormals(xa,ya);


   N = numel(xa); 
   Nt = numel(k);

   if (debug)
       plot(xa(:),ya(:),"ro");
       hold on
       plot(xa(k),ya(k),"k+");
       quiver(xa(k),ya(k),n(k,1),n(k,2),"b");
       daspect([1,1,1]);
       hold off
   end

   R0 = zeros(N);
   b0 = zeros(N,1);
   for i=1:N
       x0 = xa(i);
       y0 = ya(i);
       F = MyFun(x0,y0);
       R0(i,:) = GetRBFBase(x0,y0,xa,ya);
       b0(i) = F;
   end

   Rdb = zeros(N,Nt);
   for i=1:N
       x0 = xa(i);
       y0 = ya(i);
       [~,dRdx,dRdy] = GetRBFBase(x0,y0,xa(k),ya(k));
       Rdb(i,:) = (dRdx(:).*n(k,1)+dRdy(:).*n(k,2))';
   end

   Rc = zeros(Nt);
   bc = zeros(Nt,1);
   for i=1:Nt
       x0 = xa(k(i));
       y0 = ya(k(i));
       [~,dFdx,dFdy] = MyFun(x0,y0);
       [~,~,~,dRdxx,dRdxy,dRdyy] = GetRBFBase(x0,y0,xa(k),ya(k));
       Rc(i,:) = (dRdxx(:).*n(k,1).*n(k(i),1)+dRdxy(:).*n(k,2).*n(k(i),1)+...
                  dRdxy(:).*n(k,1).*n(k(i),2)+dRdyy(:).*n(k,2).*n(k(i),2))';
       bc(i) = dFdx*n(k(i),1)+dFdy*n(k(i),2);
   end
    
   R = [R0,Rdb;Rdb',Rc];
   b = [b0;bc];

   
   R = R0;
   b = b0;

   alpha = R\b;

   [X,Y] = meshgrid(linspace(0,1));
   Z = X*0;
   Z0 = Z;
   Res = Z;
   NN = numel(X);
   for i=1:NN
       x0 = X(i);
       y0 = Y(i);
       F = MyFun(x0,y0);
       R = GetRBFBase(x0,y0,xa,ya);
       [~,dRdx,dRdy] = GetRBFBase(x0,y0,xa(k),ya(k));
       Rbase = [R,(dRdx(:).*n(k,1)+dRdy(:).*n(k,2))'];
       Rbase = R;
       val = Rbase*alpha;
       Z(i) = val;
       Z0(i) = F;
       Res(i) = (Z0(i)-Z(i));
   end
   surf(X,Y,Z);
   fprintf("MaxAbsErr = %f\n",max(Res(:)));


end

function [F,dFdx,dFdy] = MyFun(x,y)
    F = log(1+x)+log(1+y);
    dFdx = 1/(1+x);
    dFdy = 1/(1+y);
end

function [R,dRdx,dRdy,dRdxx,dRdxy,dRdyy,LR] = GetRBFBase(x0,y0,x,y)
    dx = x-x0;
    dy = y-y0;
    r2 = dx.^2+dy.^2;
    r = sqrt(r2);
    sind = find(r<1e-8);

    R = r2(:)'.*r(:)';
    dRdx = 3*dx(:)'.*r(:)';
    dRdy = 3*dy(:)'.*r(:)';
    dRdxx = 3*(dx(:)').^2./r(:)'+3*r(:)';
    dRdxy = 3*(dx(:)'.*dy(:)')./r(:)';
    dRdyy = 3*(dy(:)').^2./r(:)'+3*r(:)';
    LR = 9*r(:)'; 

    dRdxx(sind) = 0;
    dRdxy(sind) = 0;
    dRdyy(sind) = 0;
end

function [k,n] = GetBoundaryAndNormals(x,y)
    k = boundary(x(:),y(:),0.9);
    k = k(2:end);
    N = numel(k);
    n = zeros(numel(x),2);
    for i=1:N
        ia = i-1;
        ib = i+1;
        if(ia<1)
           ia = N;
        elseif(ib>N)
           ib = 1;
        end
        ta = [x(k(ib))-x(k(ia)),y(k(ib))-y(k(ia))];
        na = [ta(2),-ta(1)]/norm(ta);
        n(k(i),:) = na;
    end
end