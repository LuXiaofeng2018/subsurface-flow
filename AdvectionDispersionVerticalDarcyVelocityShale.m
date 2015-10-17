D = .00000001;
L = .01; % in meters
C0 = 100; % concentration at base of column (Pa)
Cout = 0; % fixed percentage concentration leaving top of column (%)

rho = 2335;  % density of rock ( shale = 2335 )  (kg / m^3)
g = 9.8; % m^2 / s
eta = .05; % porosity (shale = 0.0 - 0.1) (unitless)
K = 10^-14; % hydraulic conductivity ( shale =  10-15 to 10-13 )  ( m / s)

% Pmax, pressure at the bottom of the column
% according to Wikipedia, fracking equipment can reach up to 100 megapascals (100000000 Pa) (15,000 psi)

Pmin = L * rho * g + 0; % minimum P at surface to avoid downward motion
Pmax = 10000; % some maximum P at bottom of  column

timescale = 60*60*24*365; % size of timestep , using 1 year

% simulating nonlinear change in P
%dPdX = @(x) -2 * (Pmax - Pmin) / L^2 * x;

% linear change in P
dPdX = @(x) (Pmin - Pmax) / L;  % just the slope

% calculate darcy velocity corrected for porosity
Vd = @(x) - eta * K * timescale * ( 1 + dPdX(x) / (rho * g) );  % m / timescale

%pdefun = @(x,t,u,DuDx) deal(1/D, DuDx, - ( Vd(x) / D) * DuDx);
pdefun = @(x,t,u,DuDx) deal(1, D * DuDx - Vd(x) * u, 0);  % advection and dispersion, velocity by Darcy

icfun = @(x) 0;
bcfun = @(xl,ul,xr,ur,t) deal(ul - C0 * exp(-t/10), 0, Vd(xr) * ur, 1);

xmesh = linspace(0,L,100);

years = 10;
days = 365 * years;
tspan = linspace(0,years,years * 3);  % time is in years now

m = 0;
%options = odeset('MaxStep', .0001);
%options = odeset('RelTol', 1e-8);
sol = pdepe(m,pdefun,icfun,bcfun,xmesh,tspan);
c = sol(:,:,1);
figure; surf(xmesh,tspan,c);