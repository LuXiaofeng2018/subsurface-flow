D = .3;
v = .3;
L = 20;
C0 = 100;

pdefun = @(x,t,u,DuDx) deal(1/D, DuDx, - (v / D) * DuDx);
icfun = @(x) 0;
bcfun = @(xl,ul,xr,ur,t) deal(ul - C0, 0, 0, 1);

x = linspace(0,L,30);
t = linspace(0,100,40);

m = 0;
sol = pdepe(m,pdefun,icfun,bcfun,x,t);
c = sol(:,:,1);
figure; surf(x,t,c);